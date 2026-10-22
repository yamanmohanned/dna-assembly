"""Validation of the overlap graph, the objective identities and the assembly layout (Phase 3).

Plain-language story: we made the reads ourselves, so we know which pairs really overlap. We check
that the overlap code finds every one of them (it must, for error-free reads), how many false
overlaps sneak in for different minimum-overlap settings, that the objective function and the
travelling-salesman cost agree, and that gluing the reads in their true order gives back pieces of
the real genome. Then we add typing mistakes to show why the error-tolerant overlap is needed.

Pass/fail rules are in ``configs/phase3_overlap.yaml``; results go to
``results/phase3/overlap_validation.json``.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from genomeasm import objective, provenance
from genomeasm.assembly import merged_string
from genomeasm.checks import Checks, mean_and_se
from genomeasm.checks import add_check as _add
from genomeasm.instance import AssemblyInstance
from genomeasm.io import Genome
from genomeasm.overlap import (
    expected_spurious_pairs,
    is_substring_free,
    overlap_matrix,
    overlap_matrix_naive,
    overlap_matrix_tolerant,
    truth_agreement,
)
from genomeasm.simulate import ReadLengthModel, ReadSet, SimulationConfig, simulate_reads


def _readsets(
    genome: Genome,
    coverage: float,
    length: int,
    seeds: list[np.random.SeedSequence],
    **rates: float,
) -> list[ReadSet]:
    return [
        simulate_reads(
            genome,
            SimulationConfig(
                coverage=coverage,
                read_length=ReadLengthModel.fixed(length),
                seed=int(child.generate_state(1)[0]),
                **rates,
            ),
        )
        for child in seeds
    ]


def _evaluate(
    genome: Genome,
    reads: ReadSet,
    tau: int | str,
    cost_limit: int,
    rng: np.random.Generator,
    keep_ov: bool = False,
) -> dict[str, Any]:
    """Build the instance for one read set and test every claim that can be tested on it."""
    started = time.perf_counter()
    inst = AssemblyInstance.from_readset(reads, min_overlap=tau)
    build_seconds = time.perf_counter() - started
    assert inst.true_starts is not None and inst.true_spans is not None
    truth = truth_agreement(inst.ov, inst.true_starts, inst.true_spans, inst.min_overlap)
    order = inst.true_order()
    random_order = rng.permutation(inst.n)
    asm = inst.assemble(order)
    length_true = inst.merged_length(order)
    merge_ok = len(
        merged_string(order, inst.reads, inst.ov)
    ) == length_true == asm.total_length and length_true == int(inst.lengths.sum()) - inst.fitness(
        order
    )
    cycle_ok: bool | None = None
    if inst.n <= cost_limit:
        cost = inst.cost()
        cycle_ok = all(
            objective.cycle_cost(cost, objective.tour_from_order(o)) == inst.merged_length(o)
            for o in (order, random_order)
        )
    result: dict[str, Any] = {
        "n_input": inst.n_input,
        "n": inst.n,
        "n_duplicates": inst.n_duplicates,
        "n_contained": inst.n_contained,
        "tau": inst.min_overlap,
        "build_seconds": build_seconds,
        "true_pairs": truth["true_pairs"],
        "recovered": truth["recovered"],
        "exact": truth["exact"],
        "longer": truth["longer"],
        "edges": truth["edges"],
        "spurious": truth["spurious"],
        "spurious_bound": expected_spurious_pairs(inst.n, inst.min_overlap),
        "F_true_order": inst.fitness(order),
        "L_true_order": length_true,
        "F_random_order": inst.fitness(random_order),
        "contigs_true_order": asm.n_contigs,
        "contigs_in_genome": all(s in genome.sequence for s in asm.sequences),
        "merge_identity_ok": merge_ok,
        "cycle_identity_ok": cycle_ok,
        "hist_true": truth["hist_true"],
        "hist_other": truth["hist_other"],
    }
    if keep_ov:
        result["_instance"] = inst
    return result


def _aggregate(runs: list[dict[str, Any]]) -> dict[str, Any]:
    def stat(key: str) -> dict[str, float]:
        mean, se = mean_and_se([float(r[key]) for r in runs])
        return {"mean": mean, "se": se}

    recall = [r["recovered"] / r["true_pairs"] for r in runs]
    precision = [r["exact"] / r["edges"] if r["edges"] else 1.0 for r in runs]
    recall_mean, recall_se = mean_and_se(recall)
    precision_mean, precision_se = mean_and_se(precision)
    return {
        "runs": len(runs),
        "tau": runs[0]["tau"] if len({r["tau"] for r in runs}) == 1 else [r["tau"] for r in runs],
        "n": stat("n"),
        "n_duplicates": stat("n_duplicates"),
        "n_contained": stat("n_contained"),
        "true_pairs": stat("true_pairs"),
        "recall": {"mean": recall_mean, "se": recall_se},
        "exact_share_of_edges": {"mean": precision_mean, "se": precision_se},
        "edges": stat("edges"),
        "spurious": stat("spurious"),
        "spurious_bound": stat("spurious_bound"),
        "longer_than_true": stat("longer"),
        "F_true_order": stat("F_true_order"),
        "L_true_order": stat("L_true_order"),
        "contigs_true_order": stat("contigs_true_order"),
        "contigs_in_genome_fraction": float(np.mean([r["contigs_in_genome"] for r in runs])),
        "build_seconds": stat("build_seconds"),
    }


def _main_block(
    cfg: dict[str, Any],
    genomes: dict[str, Genome],
    seeds: np.random.SeedSequence,
    checks: Checks,
    say: Callable[[str], None],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    length = int(cfg["read_length"])
    replicates = int(cfg["replicates"])
    cost_limit = int(cfg["cost_matrix_max_reads"])
    rng = np.random.default_rng(seeds.spawn(1)[0])
    rows: list[dict[str, Any]] = []
    figure_data: dict[str, Any] = {"overlap_length": [], "heatmap": None}
    heat = cfg["figure_heatmap"]
    for accession in cfg["genomes"]:
        genome = genomes[accession]
        for c in cfg["coverages"]:
            say(f"overlap graph: {accession} c={c}")
            readsets = _readsets(genome, float(c), length, seeds.spawn(replicates))
            for tau in cfg["min_overlaps"]:
                want_heat = (
                    tau == "auto" and accession == heat["genome"] and float(c) == heat["coverage"]
                )
                runs = [
                    _evaluate(genome, rs, tau, cost_limit, rng, keep_ov=want_heat and i == 0)
                    for i, rs in enumerate(readsets)
                ]
                tag = f"{accession}.c{c}.tau{tau}"
                aggregate = _aggregate(runs)
                rows.append(
                    {
                        "accession": accession,
                        "genome_length": genome.length,
                        "coverage": float(c),
                        "tau_setting": tau,
                        **aggregate,
                    }
                )
                _add(
                    checks,
                    f"overlap.recall_error_free.{tag}",
                    "every true overlap of error-free reads is found (ov >= true overlap)",
                    all(r["recovered"] == r["true_pairs"] for r in runs),
                    true_pairs=sum(r["true_pairs"] for r in runs),
                    recovered=sum(r["recovered"] for r in runs),
                )
                _add(
                    checks,
                    f"objective.identities.{tag}",
                    "L = sum(len) - F = length of the plain merge = total contig length; "
                    "cycle cost with the dummy node = L (when the matrix fits in memory)",
                    all(r["merge_identity_ok"] for r in runs)
                    and all(r["cycle_identity_ok"] is not False for r in runs),
                    cycle_identity_checked=sum(r["cycle_identity_ok"] is not None for r in runs),
                )
                if tau == "auto":
                    _add(
                        checks,
                        f"assembly.true_order_contigs_in_genome.{accession}.c{c}",
                        "gluing error-free reads in true order gives exact pieces of the genome",
                        all(r["contigs_in_genome"] for r in runs),
                    )
                    if want_heat:
                        inst = runs[0]["_instance"]
                        keep = np.flatnonzero(inst.true_starts < heat["max_start"])
                        keep = keep[np.argsort(inst.true_starts[keep], kind="stable")]
                        figure_data["heatmap"] = {
                            "accession": accession,
                            "coverage": float(c),
                            "tau": inst.min_overlap,
                            "n_shown": int(keep.size),
                            "true_starts": inst.true_starts[keep].tolist(),
                            "matrix": inst.ov[np.ix_(keep, keep)].tolist(),
                        }
                    if float(c) == float(cfg["figure_overlap_length"]["coverage"]):
                        figure_data["overlap_length"].append(
                            {
                                "accession": accession,
                                "coverage": float(c),
                                "tau": runs[0]["tau"],
                                "runs": replicates,
                                "hist_true": np.sum(
                                    [r["hist_true"] for r in runs], axis=0
                                ).tolist(),
                                "hist_other": np.sum(
                                    [r["hist_other"] for r in runs], axis=0
                                ).tolist(),
                            }
                        )
                for r in runs:
                    r.pop("_instance", None)
    return rows, figure_data


def _equivalence_block(
    cfg: dict[str, Any], genomes: dict[str, Genome], seeds: np.random.SeedSequence, checks: Checks
) -> dict[str, Any]:
    block = cfg["equivalence"]
    genome = genomes[block["genome"]]
    rs = _readsets(genome, float(block["coverage"]), int(cfg["read_length"]), seeds.spawn(1))[0]
    inst = AssemblyInstance.from_readset(rs, min_overlap=int(block["min_overlap"]))
    tau = int(block["min_overlap"])
    naive = overlap_matrix_naive(inst.reads, tau)
    fast = overlap_matrix(inst.reads, tau)
    tolerant_zero = overlap_matrix_tolerant(inst.reads, tau, 0.0)
    _add(
        checks,
        "overlap.fast_equals_naive",
        "prefix-index overlap code equals the naive definition on every pair",
        np.array_equal(naive, fast),
        n_reads=inst.n,
        n_pairs=inst.n * (inst.n - 1),
    )
    _add(
        checks,
        "overlap.tolerant_with_zero_tolerance_equals_exact",
        "error-tolerant code with tolerance 0 equals the exact overlap code",
        np.array_equal(tolerant_zero, fast),
    )
    _add(
        checks,
        "overlap.preprocessed_reads_are_substring_free",
        "no kept read is identical to or inside another kept read",
        is_substring_free(inst.reads),
        n_reads=inst.n,
        removed_duplicates=inst.n_duplicates,
        removed_contained=inst.n_contained,
    )
    return {"n_reads": inst.n, "min_overlap": tau}


def _noise_block(
    cfg: dict[str, Any], genomes: dict[str, Genome], seeds: np.random.SeedSequence, checks: Checks
) -> dict[str, Any]:
    block = cfg["noise"]
    tol = cfg["tolerances"]
    genome = genomes[block["genome"]]
    readsets = _readsets(
        genome,
        float(block["coverage"]),
        int(cfg["read_length"]),
        seeds.spawn(int(block["replicates"])),
        substitution_rate=float(block["substitution_rate"]),
    )
    out: dict[str, Any] = {"config": dict(block), "modes": {}}
    for mode, tolerance in (("exact", 0.0), ("tolerant", float(block["tolerance"]))):
        recalls: list[float] = []
        edges: list[float] = []
        spurious: list[float] = []
        for rs in readsets:
            inst = AssemblyInstance.from_sequences(
                rs.sequences,
                min_overlap=int(block["min_overlap"]),
                tolerance=tolerance,
                starts=rs.starts,
                spans=rs.spans,
            )
            assert inst.true_starts is not None and inst.true_spans is not None
            truth = truth_agreement(inst.ov, inst.true_starts, inst.true_spans, inst.min_overlap)
            recalls.append(truth["recovered"] / truth["true_pairs"])
            edges.append(truth["edges"])
            spurious.append(truth["spurious"])
        recall_mean, recall_se = mean_and_se(recalls)
        out["modes"][mode] = {
            "tolerance": tolerance,
            "recall": {"mean": recall_mean, "se": recall_se},
            "edges": mean_and_se(edges)[0],
            "spurious_edges": mean_and_se(spurious)[0],
        }
    exact_recall = out["modes"]["exact"]["recall"]["mean"]
    tolerant_recall = out["modes"]["tolerant"]["recall"]["mean"]
    _add(
        checks,
        "noise.exact_overlaps_lose_true_pairs",
        "with 1% substitution errors exact overlaps miss many true overlaps",
        exact_recall <= tol["noise_exact_recall_max"],
        observed=exact_recall,
        limit=tol["noise_exact_recall_max"],
    )
    _add(
        checks,
        "noise.tolerant_overlaps_recover_true_pairs",
        "error-tolerant overlaps recover most true overlaps under the same errors",
        tolerant_recall >= tol["noise_tolerant_recall_min"],
        observed=tolerant_recall,
        limit=tol["noise_tolerant_recall_min"],
    )
    return out


def run_overlap_validation(
    cfg: dict[str, Any],
    genomes: dict[str, Genome],
    *,
    repo_dir: Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run all Phase 3 checks and return the result dictionary."""
    say = progress or (lambda message: None)
    root = np.random.SeedSequence(int(cfg["seed"]))
    s_main, s_equiv, s_noise = root.spawn(3)
    checks: Checks = []
    rows, figure_data = _main_block(cfg, genomes, s_main, checks, say)
    say("equivalence of fast and naive code")
    equivalence = _equivalence_block(cfg, genomes, s_equiv, checks)
    say("noise and error-tolerant overlaps")
    noise = _noise_block(cfg, genomes, s_noise, checks)
    return {
        "name": cfg["name"],
        "config": cfg,
        "provenance": provenance.collect(cfg, seed=int(cfg["seed"]), repo_dir=repo_dir),
        "instances": rows,
        "equivalence": equivalence,
        "noise": noise,
        "figure_data": figure_data,
        "checks": checks,
        "n_checks": len(checks),
        "n_failed": sum(not c["passed"] for c in checks),
        "all_passed": all(c["passed"] for c in checks),
    }
