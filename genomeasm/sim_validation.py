"""Statistical validation of the simulator against theory (Phase 2 exit gate).

Plain-language story: before we trust the simulated machine we test it like a measuring
instrument. We run it many times (Monte Carlo) and compare what it produces with what
mathematics says it must produce: how many positions stay uncovered, how many separate
contigs appear, how often it makes typing mistakes, and whether read lengths and start
positions follow their rules. Each comparison is a *check* with a tolerance that is written
in ``configs/phase2_validation.yaml`` before the run.

All randomness comes from ``numpy.random.SeedSequence(seed)`` children, so a run is reproducible.
A check passes when the observed value is within ``z_max`` standard errors of the expected value
(or the other documented rule). Standard error = standard deviation across replicates divided by
the square root of the number of replicates.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

from genomeasm import coverage as cov
from genomeasm import provenance
from genomeasm.checks import Checks, mean_and_se
from genomeasm.checks import add_check as _add
from genomeasm.io import Genome
from genomeasm.simulate import (
    ReadLengthModel,
    SimulationConfig,
    sample_geometry,
    simulate_reads,
)


def edit_distance(a: str, b: str) -> int:
    """Levenshtein distance (insert, delete, substitute each cost 1), classic dynamic programme."""
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def _within(observed: float, expected: float, se: float, z_max: float, floor: float = 0.0) -> bool:
    return abs(observed - expected) <= max(z_max * se, floor)


def histogram_z_scores(
    histograms: np.ndarray,
    exact: np.ndarray,
    genome_length: int,
    min_positions: float,
) -> tuple[list[tuple[int, int]], np.ndarray]:
    """z-score of the simulated depth histogram against the exact one, per group of bins.

    Neighbouring bins are merged until each group expects at least ``min_positions`` positions
    (``coverage.group_bins``), because a bin that never occurs in any run has no measurable
    spread. For each group, z = (mean over runs - exact) / (std over runs / sqrt(runs)).
    ``histograms`` has shape (runs, bins).

    History (decision D-022): the first version of this check compared every single bin with a
    tiny absolute floor and failed on sparse bins that never occurred in the runs (zero spread).
    A 2000-run diagnostic showed no real discrepancy, so the method was changed, not the tolerance.
    """
    groups = cov.group_bins(exact, genome_length, min_positions)
    observed = cov.grouped_fractions(histograms, groups)
    expected = cov.grouped_fractions(exact, groups)
    runs = observed.shape[0]
    diff = observed.mean(axis=0) - expected
    se = observed.std(axis=0, ddof=1) / np.sqrt(runs)
    # A spread below 1e-12 is floating-point noise, not data (fractions are multiples of 1/G).
    measurable = se > 1e-12
    safe = np.where(measurable, se, 1.0)
    z = np.where(measurable, diff / safe, np.where(np.abs(diff) < 1e-9, 0.0, np.inf))
    return groups, z


# --------------------------------------------------------------------------------------
# Coverage and contigs on the real genomes
# --------------------------------------------------------------------------------------


def _coverage_block(
    cfg: dict[str, Any],
    genomes: dict[str, Genome],
    seeds: np.random.SeedSequence,
    checks: Checks,
) -> list[dict[str, Any]]:
    tol = cfg["tolerances"]
    length = int(cfg["read_length"])
    theta = int(cfg["min_overlap"])
    n_bins = int(cfg["depth_bins"])
    replicates = int(cfg["replicates"])
    out: list[dict[str, Any]] = []
    for accession in cfg["genomes"]:
        genome = genomes[accession]
        g = genome.length
        for c in cfg["coverages"]:
            sim = SimulationConfig(coverage=float(c), read_length=ReadLengthModel.fixed(length))
            n_reads = sim.n_reads(g)
            c_nominal = sim.nominal_coverage(g)
            support, pmf = sim.read_length.pmf()
            uncovered = np.empty(replicates)
            contigs = np.empty(replicates)
            hists = np.empty((replicates, n_bins))
            identity_ok = True
            for r, child in enumerate(seeds.spawn(replicates)):
                starts, spans = sample_geometry(g, sim, np.random.default_rng(child))
                depth = cov.coverage_depth(g, starts, spans)
                identity_ok &= int(depth.sum()) == int(spans.sum()) == n_reads * length
                uncovered[r] = cov.uncovered_fraction(depth)
                contigs[r] = cov.count_contigs(starts, spans, theta)
                hists[r] = cov.depth_histogram(depth, n_bins)
            exact_u = cov.expected_uncovered_fraction(g, n_reads, support, pmf)
            lw_u = cov.lw_uncovered_fraction(c_nominal)
            exact_h = cov.expected_depth_histogram(g, n_reads, support, pmf, n_bins)
            lw_h = cov.lw_depth_pmf(c_nominal, n_bins)
            exact_k = cov.expected_contigs_exact(g - length + 1, n_reads, length - theta)
            lw_k = cov.lw_expected_contigs(n_reads, c_nominal, theta, length)
            u_mean, u_se = mean_and_se(uncovered)
            k_mean, k_se = mean_and_se(contigs)
            h_mean = hists.mean(axis=0)
            tag = f"{accession}.c{c}"

            _add(
                checks,
                f"coverage.identity.{tag}",
                "sum of depths equals N*l exactly, so realised coverage equals N*l/G",
                identity_ok,
            )
            _add(
                checks,
                f"coverage.uncovered_vs_exact.{tag}",
                "mean uncovered fraction agrees with the exact finite-genome formula",
                _within(u_mean, exact_u, u_se, tol["z_max"]),
                observed=u_mean,
                expected=exact_u,
                standard_error=u_se,
                z_max=tol["z_max"],
            )
            groups, z_hist = histogram_z_scores(hists, exact_h, g, tol["histogram_min_positions"])
            max_z = float(np.max(np.abs(z_hist)))
            _add(
                checks,
                f"coverage.histogram_vs_exact.{tag}",
                "depth histogram agrees with the exact Binomial mixture in every bin group",
                max_z <= tol["histogram_z_max"],
                max_abs_z=max_z,
                n_groups=len(groups),
                histogram_z_max=tol["histogram_z_max"],
                max_abs_difference=float(np.max(np.abs(h_mean - exact_h))),
            )
            _add(
                checks,
                f"contigs.vs_exact.{tag}",
                "mean number of true contigs agrees with the derived exact formula",
                _within(k_mean, exact_k, k_se, tol["z_max"], floor=tol["contig_floor"]),
                observed=k_mean,
                expected=exact_k,
                standard_error=k_se,
                z_max=tol["z_max"],
            )
            out.append(
                {
                    "accession": accession,
                    "genome_length": g,
                    "coverage_requested": float(c),
                    "n_reads": n_reads,
                    "read_length": length,
                    "nominal_coverage": c_nominal,
                    "replicates": replicates,
                    "min_overlap": theta,
                    "uncovered": {
                        "mean": u_mean,
                        "se": u_se,
                        "q25": float(np.quantile(uncovered, 0.25)),
                        "q75": float(np.quantile(uncovered, 0.75)),
                        "exact": exact_u,
                        "lander_waterman": lw_u,
                    },
                    "contigs": {
                        "mean": k_mean,
                        "se": k_se,
                        "exact": exact_k,
                        "lander_waterman": lw_k,
                    },
                    "histogram": {
                        "bins": n_bins,
                        "observed_mean": h_mean.tolist(),
                        "observed_q25": np.quantile(hists, 0.25, axis=0).tolist(),
                        "observed_q75": np.quantile(hists, 0.75, axis=0).tolist(),
                        "exact": exact_h.tolist(),
                        "poisson": lw_h.tolist(),
                    },
                }
            )
    return out


def _large_genome_block(
    cfg: dict[str, Any], seeds: np.random.SeedSequence, checks: Checks
) -> list[dict[str, Any]]:
    """Lander-Waterman e^-c is an approximation for a long genome: test it where it applies."""
    tol = cfg["tolerances"]
    block = cfg["large_genome_check"]
    g = int(block["genome_length"])
    length = int(cfg["read_length"])
    replicates = int(block["replicates"])
    out: list[dict[str, Any]] = []
    for c in block["coverages"]:
        sim = SimulationConfig(coverage=float(c), read_length=ReadLengthModel.fixed(length))
        c_nominal = sim.nominal_coverage(g)
        values = np.empty(replicates)
        for r, child in enumerate(seeds.spawn(replicates)):
            starts, spans = sample_geometry(g, sim, np.random.default_rng(child))
            values[r] = cov.uncovered_fraction(cov.coverage_depth(g, starts, spans))
        mean, se = mean_and_se(values)
        lw = cov.lw_uncovered_fraction(c_nominal)
        allowed = tol["z_max"] * se + tol["lw_large_genome_rel"] * lw
        _add(
            checks,
            f"coverage.uncovered_vs_lw_large_genome.c{c}",
            f"on a {g:,}-position genome the uncovered fraction matches exp(-c)",
            abs(mean - lw) <= allowed,
            observed=mean,
            expected=lw,
            standard_error=se,
            allowed_difference=allowed,
        )
        out.append(
            {
                "genome_length": g,
                "coverage_requested": float(c),
                "nominal_coverage": c_nominal,
                "replicates": replicates,
                "uncovered_mean": mean,
                "uncovered_se": se,
                "lander_waterman": lw,
                "relative_gap": mean / lw - 1.0,
            }
        )
    return out


# --------------------------------------------------------------------------------------
# Error model
# --------------------------------------------------------------------------------------


def _rate_check(
    checks: Checks,
    check_id: str,
    description: str,
    events: int,
    trials: int,
    p: float,
    z_max: float,
) -> dict[str, float]:
    observed = events / trials
    se = float(np.sqrt(p * (1.0 - p) / trials))
    _add(
        checks,
        check_id,
        description,
        abs(observed - p) <= z_max * se,
        observed=observed,
        expected=p,
        standard_error=se,
        z_max=z_max,
        trials=trials,
    )
    return {"observed": observed, "expected": p, "standard_error": se, "trials": trials}


def _error_block(
    cfg: dict[str, Any],
    genomes: dict[str, Genome],
    seeds: np.random.SeedSequence,
    checks: Checks,
) -> dict[str, Any]:
    tol = cfg["tolerances"]
    block = cfg["error_check"]
    genome = genomes[block["genome"]]
    replicates = int(block["replicates"])
    length = int(cfg["read_length"])
    out: dict[str, Any] = {"substitution_only": [], "indel": {}}

    def one_run(child: np.random.SeedSequence, **rates: float) -> Any:
        seed = int(child.generate_state(1)[0])
        sim = SimulationConfig(
            coverage=float(block["coverage"]),
            read_length=ReadLengthModel.fixed(length),
            seed=seed,
            **rates,
        )
        return simulate_reads(genome, sim)

    for p in cfg["substitution_rates"]:
        z_scores: list[float] = []
        total_sub = total_bases = 0
        hamming_matches_counts = True
        for child in seeds.spawn(replicates):
            rs = one_run(child, substitution_rate=float(p))
            bases = int(rs.spans.sum())
            subs = int(rs.n_substitutions.sum())
            total_sub += subs
            total_bases += bases
            z_scores.append((subs - p * bases) / np.sqrt(p * (1 - p) * bases))
            for i in range(0, rs.n_reads, 7):  # independent recount on every 7th read
                window = genome.sequence[rs.starts[i] : rs.starts[i] + rs.spans[i]]
                hamming = sum(x != y for x, y in zip(rs.sequences[i], window, strict=True))
                hamming_matches_counts &= hamming == int(rs.n_substitutions[i])
        summary = _rate_check(
            checks,
            f"errors.substitution_rate.p{p}",
            "pooled substitution rate equals the configured rate",
            total_sub,
            total_bases,
            float(p),
            tol["z_max"],
        )
        z = np.asarray(z_scores)
        z_mean_limit = tol["z_max"] / np.sqrt(z.size)
        lo, hi = tol["z_variance_range"]
        _add(
            checks,
            f"errors.substitution_calibration.p{p}",
            "per-run z-scores have mean ~0 and variance ~1 (unbiased rate, correct noise)",
            abs(z.mean()) <= z_mean_limit and lo <= z.var(ddof=1) <= hi,
            z_mean=z.mean(),
            z_variance=z.var(ddof=1),
            z_mean_limit=z_mean_limit,
            z_variance_range=[lo, hi],
        )
        _add(
            checks,
            f"errors.hamming_equals_recorded.p{p}",
            "substitution-only reads: Hamming distance to the true window equals recorded count",
            hamming_matches_counts,
        )
        out["substitution_only"].append({"rate": float(p), **summary, "runs": replicates})

    scenario = cfg["indel_scenario"]
    p_sub, p_ins, p_del = (
        scenario[k] for k in ("substitution_rate", "insertion_rate", "deletion_rate")
    )
    spans_total = subs_total = ins_total = del_total = 0
    lengths_ok = True
    edit_bound_ok = True
    checked_edit = 0
    for run, child in enumerate(seeds.spawn(replicates)):
        rs = one_run(child, **scenario)
        spans_total += int(rs.spans.sum())
        subs_total += int(rs.n_substitutions.sum())
        ins_total += int(rs.n_insertions.sum())
        del_total += int(rs.n_deletions.sum())
        lengths_ok &= bool(np.all(rs.lengths == rs.spans - rs.n_deletions + rs.n_insertions))
        if run < 3:  # edit distance is O(l^2) per read, so only a sample is checked
            for i in range(0, rs.n_reads, 11):
                window = genome.sequence[rs.starts[i] : rs.starts[i] + rs.spans[i]]
                edits = int(rs.n_substitutions[i] + rs.n_insertions[i] + rs.n_deletions[i])
                edit_bound_ok &= edit_distance(rs.sequences[i], window) <= edits
                checked_edit += 1
    out["indel"] = {
        "scenario": dict(scenario),
        "deletion": _rate_check(
            checks,
            "errors.deletion_rate",
            "pooled deletion rate equals the configured rate",
            del_total,
            spans_total,
            float(p_del),
            tol["z_max"],
        ),
        "insertion": _rate_check(
            checks,
            "errors.insertion_rate",
            "pooled insertion rate equals the configured rate",
            ins_total,
            spans_total,
            float(p_ins),
            tol["z_max"],
        ),
        "substitution_among_kept": _rate_check(
            checks,
            "errors.substitution_rate_with_indels",
            "substitution rate among kept bases equals the configured rate",
            subs_total,
            spans_total - del_total,
            float(p_sub),
            tol["z_max"],
        ),
        "reads_edit_checked": checked_edit,
    }
    _add(
        checks,
        "errors.read_length_identity",
        "read length equals span - deletions + insertions for every read",
        lengths_ok,
    )
    _add(
        checks,
        "errors.edit_distance_bound",
        "edit distance between read and true window never exceeds the recorded number of edits",
        edit_bound_ok,
        reads_checked=checked_edit,
    )
    return out


# --------------------------------------------------------------------------------------
# Read lengths and start positions
# --------------------------------------------------------------------------------------


def _length_block(
    cfg: dict[str, Any], seeds: np.random.SeedSequence, checks: Checks
) -> list[dict[str, Any]]:
    tol = cfg["tolerances"]
    block = cfg["length_check"]
    g = int(block["genome_length"])
    out: list[dict[str, Any]] = []
    models: list[tuple[str, ReadLengthModel]] = [
        ("fixed", ReadLengthModel.fixed(int(cfg["read_length"])))
    ]
    for spec in cfg["length_models"]:
        spec = dict(spec)
        name = spec.pop("name")
        models.append((name, ReadLengthModel.from_dict(spec)))
    for name, model in models:
        coverage = block["n_reads"] * model.mean / g
        sim = SimulationConfig(coverage=coverage, read_length=model)
        rng = np.random.default_rng(seeds.spawn(1)[0])
        starts, spans = sample_geometry(g, sim, rng)
        n = starts.size
        support, pmf = model.pmf()
        mean_expected = model.mean
        sd_expected = float(np.sqrt(np.dot(pmf, (support - mean_expected) ** 2)))
        mean_obs = float(spans.mean())
        se = sd_expected / np.sqrt(n)
        in_support = bool(np.all((spans >= model.low) & (spans <= model.high)))
        counts = np.bincount(spans - model.low, minlength=support.size).astype(np.float64)
        if support.size > 1:
            chi = stats.chisquare(counts, pmf * n)
            p_value = float(chi.pvalue)
            gof_ok = p_value >= tol["chi2_alpha"]
        else:
            p_value, gof_ok = 1.0, bool(counts[0] == n)
        u = starts / (g - spans + 1)  # position in [0,1) among the allowed start positions
        ks = stats.kstest(u, "uniform")
        _add(
            checks,
            f"lengths.support.{name}",
            "all read lengths lie in the model's range",
            in_support,
        )
        _add(
            checks,
            f"lengths.mean.{name}",
            "mean read length agrees with the model's exact mean",
            _within(mean_obs, mean_expected, se, tol["z_max"]),
            observed=mean_obs,
            expected=mean_expected,
            standard_error=se,
        )
        _add(
            checks,
            f"lengths.distribution.{name}",
            "read-length frequencies follow the model (chi-square goodness of fit)",
            gof_ok,
            p_value=p_value,
            alpha=tol["chi2_alpha"],
        )
        _add(
            checks,
            f"starts.uniform.{name}",
            "start positions are uniform among the positions where the read fits (KS test)",
            float(ks.pvalue) >= tol["chi2_alpha"],
            p_value=float(ks.pvalue),
            alpha=tol["chi2_alpha"],
        )
        out.append(
            {
                "name": name,
                "model": model.to_dict(),
                "n_reads": int(n),
                "mean_observed": mean_obs,
                "mean_expected": mean_expected,
                "chi2_p_value": p_value,
                "ks_p_value": float(ks.pvalue),
            }
        )
    return out


# --------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------


def run_validation(
    cfg: dict[str, Any],
    genomes: dict[str, Genome],
    *,
    repo_dir: Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run every check of the Phase 2 validation and return the full result dictionary."""
    say = progress or (lambda message: None)
    root = np.random.SeedSequence(int(cfg["seed"]))
    s_cov, s_large, s_err, s_len = root.spawn(4)
    checks: Checks = []
    say("coverage and contigs on the real genomes")
    coverage = _coverage_block(cfg, genomes, s_cov, checks)
    say("Lander-Waterman on a large genome")
    large = _large_genome_block(cfg, s_large, checks)
    say("error model")
    errors = _error_block(cfg, genomes, s_err, checks)
    say("read lengths and start positions")
    lengths = _length_block(cfg, s_len, checks)
    return {
        "name": cfg["name"],
        "config": cfg,
        "provenance": provenance.collect(cfg, seed=int(cfg["seed"]), repo_dir=repo_dir),
        "coverage": coverage,
        "large_genome": large,
        "errors": errors,
        "lengths": lengths,
        "checks": checks,
        "n_checks": len(checks),
        "n_failed": sum(not c["passed"] for c in checks),
        "all_passed": all(c["passed"] for c in checks),
    }
