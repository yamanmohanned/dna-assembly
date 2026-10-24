"""Command-line interface: ``genomeasm data fetch | verify | cards``."""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Sequence
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from genomeasm import __version__, provenance
from genomeasm.datacard import write_cards
from genomeasm.instance import AssemblyInstance
from genomeasm.io import (
    GenomeError,
    GenomeNotCachedError,
    download_genome,
    fetch_from_ncbi,
    load_genome,
)
from genomeasm.overlap_validation import run_overlap_validation
from genomeasm.plots import make_phase2_figures, make_phase3_figures
from genomeasm.sim_validation import run_validation
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads
from genomeasm.solver_validation import run_solver_validation
from genomeasm.solvers.alns import alns
from genomeasm.solvers.common import Budget, Problem
from genomeasm.solvers.greedy import greedy, nearest_neighbour
from genomeasm.solvers.lkh import LkhError, lkh
from genomeasm.solvers.tabu import tabu_search
from genomeasm.solvers.warmup import warmup

DEFAULT_CONFIG = Path("configs/genomes.yaml")
DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_READS_DIR = Path("data/processed/reads")
DEFAULT_ASSEMBLY_DIR = Path("data/processed/assemblies")


def load_accessions(config_path: Path) -> list[str]:
    """Read the list of accessions from a YAML config such as configs/genomes.yaml."""
    data = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    entries = data["genomes"]
    return [str(entry["accession"]) for entry in entries]


def _accessions(args: argparse.Namespace) -> list[str]:
    return list(args.accession) if args.accession else load_accessions(args.config)


def _cmd_fetch(args: argparse.Namespace) -> int:
    email = args.email or os.environ.get("NCBI_EMAIL")
    fetcher = partial(fetch_from_ncbi, email=email)
    failures = 0
    for accession in _accessions(args):
        try:
            if not args.force:
                try:
                    loaded = load_genome(accession, args.raw_dir)
                    print(f"{accession}: already cached and verified ({loaded.genome.length} bp)")
                    continue
                except GenomeNotCachedError:
                    pass
            loaded = download_genome(accession, args.raw_dir, fetcher=fetcher)
            print(f"{accession}: downloaded and verified ({loaded.genome.length} bp)")
        except GenomeError as exc:
            failures += 1
            print(f"error: {exc}", file=sys.stderr)
    return 1 if failures else 0


def _cmd_verify(args: argparse.Namespace) -> int:
    failures = 0
    for accession in _accessions(args):
        try:
            loaded = load_genome(accession, args.raw_dir)
            print(f"{accession}: OK, {loaded.genome.length} bp, sha256 {loaded.file_sha256}")
        except GenomeError as exc:
            failures += 1
            print(f"error: {exc}", file=sys.stderr)
    return 1 if failures else 0


def _cmd_cards(args: argparse.Namespace) -> int:
    try:
        paths = write_cards(args.raw_dir, args.processed_dir, _accessions(args))
    except GenomeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for path in paths:
        print(f"wrote {path}")
    return 0


def _json_default(value: object) -> object:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialise {type(value).__name__}")


def _write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(data, indent=2, default=_json_default) + "\n"
    path.write_text(text, encoding="utf-8", newline="\n")


def _cmd_sim_run(args: argparse.Namespace) -> int:
    try:
        genome = load_genome(args.accession, args.raw_dir).genome
        if args.read_length_sd > 0:
            spread = round(3 * args.read_length_sd)
            model = ReadLengthModel.normal(
                args.read_length,
                args.read_length_sd,
                low=max(1, args.read_length - spread),
                high=args.read_length + spread,
            )
        else:
            model = ReadLengthModel.fixed(int(args.read_length))
        config = SimulationConfig(
            coverage=args.coverage,
            read_length=model,
            substitution_rate=args.sub_rate,
            insertion_rate=args.ins_rate,
            deletion_rate=args.del_rate,
            seed=args.seed,
        )
        reads = simulate_reads(genome, config)
    except (GenomeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    out_dir = (
        args.out_dir or DEFAULT_READS_DIR / f"{args.accession}_c{args.coverage:g}_s{args.seed}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "reads.fasta").write_text(reads.to_fasta(), encoding="utf-8", newline="\n")
    summary = {
        "accession": args.accession,
        "genome_length": genome.length,
        "simulation": config.to_dict(),
        "n_reads": reads.n_reads,
        "nominal_coverage": config.nominal_coverage(genome.length),
        "realised_coverage": reads.realised_coverage,
        "total_substitutions": int(reads.n_substitutions.sum()),
        "total_insertions": int(reads.n_insertions.sum()),
        "total_deletions": int(reads.n_deletions.sum()),
        "provenance": provenance.collect(config.to_dict(), seed=args.seed, repo_dir=Path.cwd()),
    }
    _write_json(out_dir / "reads.config.json", summary)
    print(
        f"{reads.n_reads} reads (coverage {summary['nominal_coverage']:.2f}x) written to {out_dir}"
    )
    return 0


def _cmd_sim_validate(args: argparse.Namespace) -> int:
    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    wanted = list(dict.fromkeys([*cfg["genomes"], cfg["error_check"]["genome"]]))
    try:
        genomes = {acc: load_genome(acc, args.raw_dir).genome for acc in wanted}
    except GenomeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    results = run_validation(cfg, genomes, repo_dir=Path.cwd(), progress=print)
    results_path = Path(args.results_dir) / "validation.json"
    _write_json(results_path, results)
    figures = make_phase2_figures(results, Path(args.figures_dir), results_path.as_posix())
    print(f"{results['n_checks']} checks, {results['n_failed']} failed; wrote {results_path}")
    print(f"wrote {len(figures)} figure files to {args.figures_dir}")
    return 0 if results["all_passed"] else 1


def _cmd_sim_figures(args: argparse.Namespace) -> int:
    results_path = Path(args.results_dir) / "validation.json"
    results = json.loads(results_path.read_text(encoding="utf-8"))
    figures = make_phase2_figures(results, Path(args.figures_dir), results_path.as_posix())
    print(f"wrote {len(figures)} figure files to {args.figures_dir}")
    return 0


def _cmd_overlap_validate(args: argparse.Namespace) -> int:
    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    wanted = list(
        dict.fromkeys(
            [*cfg["genomes"], cfg["equivalence"]["genome"], cfg["noise"]["genome"]],
        )
    )
    try:
        genomes = {acc: load_genome(acc, args.raw_dir).genome for acc in wanted}
    except GenomeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    results = run_overlap_validation(cfg, genomes, repo_dir=Path.cwd(), progress=print)
    results_path = Path(args.results_dir) / "overlap_validation.json"
    _write_json(results_path, results)
    figures = make_phase3_figures(results, Path(args.figures_dir), results_path.as_posix())
    print(f"{results['n_checks']} checks, {results['n_failed']} failed; wrote {results_path}")
    print(f"wrote {len(figures)} figure files to {args.figures_dir}")
    return 0 if results["all_passed"] else 1


def _cmd_overlap_figures(args: argparse.Namespace) -> int:
    results_path = Path(args.results_dir) / "overlap_validation.json"
    results = json.loads(results_path.read_text(encoding="utf-8"))
    figures = make_phase3_figures(results, Path(args.figures_dir), results_path.as_posix())
    print(f"wrote {len(figures)} figure files to {args.figures_dir}")
    return 0


def _cmd_assemble(args: argparse.Namespace) -> int:
    """Walking skeleton: simulate reads, build the overlap graph, order, layout, consensus."""
    try:
        genome = load_genome(args.accession, args.raw_dir).genome
        config = SimulationConfig(
            coverage=args.coverage,
            read_length=ReadLengthModel.fixed(int(args.read_length)),
            substitution_rate=args.sub_rate,
            seed=args.seed,
        )
        reads = simulate_reads(genome, config)
        tau: int | str = "auto" if args.min_overlap == "auto" else int(args.min_overlap)
        inst = AssemblyInstance.from_readset(reads, min_overlap=tau, tolerance=args.tolerance)
    except (GenomeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    if args.order == "true":
        order = inst.true_order()
    elif args.order == "random":
        order = np.random.default_rng(args.seed).permutation(inst.n).astype(np.int64)
    else:  # trivial: the reads in the order the sequencer produced them (random positions)
        order = np.arange(inst.n, dtype=np.int64)
    assembly = inst.assemble(order)
    out_dir = args.out_dir or (
        DEFAULT_ASSEMBLY_DIR / f"{args.accession}_c{args.coverage:g}_s{args.seed}_{args.order}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = "".join(
        f">contig_{i + 1:03d} length={len(c.sequence)} reads={len(c.reads)}\n{c.sequence}\n"
        for i, c in enumerate(assembly.contigs)
    )
    (out_dir / "contigs.fasta").write_text(fasta, encoding="utf-8", newline="\n")
    summary = {
        "accession": args.accession,
        "genome_length": genome.length,
        "simulation": config.to_dict(),
        "instance": inst.summary(),
        "order": args.order,
        "fitness_F": inst.fitness(order),
        "merged_length_L": inst.merged_length(order),
        "n_contigs": assembly.n_contigs,
        "contig_lengths": [len(c.sequence) for c in assembly.contigs],
        "contigs_exact_in_genome": [c.sequence in genome.sequence for c in assembly.contigs],
        "provenance": provenance.collect(
            {"simulation": config.to_dict(), "min_overlap": str(args.min_overlap),
             "tolerance": args.tolerance, "order": args.order},
            seed=args.seed,
            repo_dir=Path.cwd(),
        ),
    }  # fmt: skip
    _write_json(out_dir / "summary.json", summary)
    print(
        f"{inst.n} reads (from {inst.n_input}), tau={inst.min_overlap}, order={args.order}: "
        f"F={summary['fitness_F']}, L={summary['merged_length_L']}, "
        f"{assembly.n_contigs} contig(s), written to {out_dir}"
    )
    return 0


def _cmd_solvers_validate(args: argparse.Namespace) -> int:
    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    try:
        genomes = {cfg["genome"]: load_genome(cfg["genome"], args.raw_dir).genome}
    except GenomeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    results = run_solver_validation(cfg, genomes, repo_dir=Path.cwd(), progress=print)
    path = Path(args.results_dir) / "correctness.json"
    _write_json(path, results)
    print(f"{results['n_checks']} checks, {results['n_failed']} failed; wrote {path}")
    return 0 if results["all_passed"] else 1


def _cmd_app(args: argparse.Namespace) -> int:
    """app run: start the Streamlit application; app benchmark: time its worst interactive runs."""
    import subprocess

    from genomeasm.appkit.paths import ROOT

    if args.command == "run":
        command = [
            sys.executable, "-m", "streamlit", "run", str(ROOT / "app" / "streamlit_app.py"),
            "--server.port", str(args.port), "--server.headless", str(not args.browser).lower(),
            "--browser.gatherUsageStats", "false",
        ]  # fmt: skip
        return subprocess.call(command, cwd=ROOT)
    from genomeasm.appkit import timing

    result = timing.measure()
    result["provenance"] = provenance.collect(
        {"limits": result["interface_limits"], "scenarios": sorted(result["scenarios"])},
        seed=None, repo_dir=Path.cwd(),
    )  # fmt: skip
    _write_json(args.out, result)
    for name, scenario in result["scenarios"].items():
        print(f"{name}: {scenario['total']:.1f} s  {scenario}")
    print(f"worst case {result['worst_total_seconds']:.1f} s of {result['limit_seconds']:.0f} s")
    return 0 if result["ok"] else 1


def _cmd_docs(args: argparse.Namespace) -> int:
    """docs build: the PDF report, the slides (pptx and PDF) and docs/ORAL_EXAM.md from results/."""
    if args.command == "verify":
        from genomeasm.docgen import verify
        from genomeasm.docgen.facts import Facts

        report = verify.verify_all(Facts.load())
        _write_json(args.out, report)
        for name, doc in report["documents"].items():
            print(
                f"{name}: {doc['numbers_checked']} numbers, unexplained: {doc['unexplained']}"
            )
        print("recomputation problems:", report["recomputation_problems"])
        print("figure problems:", report["figure_problems"])
        print("OK" if report["ok"] else "NOT OK")
        return 0 if report["ok"] else 1
    from genomeasm.docgen.build import build_all

    for target, files in build_all(args.only, author=args.author).items():
        for file in files:
            print(f"{target}: {file} ({file.stat().st_size:,} bytes)")
    return 0


def _cmd_solve(args: argparse.Namespace) -> int:
    """Simulate reads, build the instance and run one solver (a quick look, not an experiment)."""
    try:
        genome = load_genome(args.accession, args.raw_dir).genome
        reads = simulate_reads(
            genome,
            SimulationConfig(
                coverage=args.coverage,
                read_length=ReadLengthModel.fixed(int(args.read_length)),
                substitution_rate=args.sub_rate,
                seed=args.read_seed,
            ),
        )
        tau: int | str = "auto" if args.min_overlap == "auto" else int(args.min_overlap)
        inst = AssemblyInstance.from_readset(reads, min_overlap=tau, tolerance=args.tolerance)
        problem = Problem.from_instance(inst)
        warmup()
        budget = Budget(seconds=args.seconds)
        if args.solver == "greedy":
            result = greedy(problem)
        elif args.solver == "nn":
            result = nearest_neighbour(problem)
        elif args.solver == "tabu":
            result = tabu_search(problem, budget, seed=args.seed)
        elif args.solver == "alns":
            result = alns(problem, budget, seed=args.seed, policy=args.policy)
        else:
            result = lkh(problem, seconds=args.seconds, seed=args.seed)
    except (GenomeError, ValueError, LkhError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(
        f"{result.solver}: n={inst.n}, F={result.fitness}, L={result.merged_length}, "
        f"F(true order)={inst.fitness(inst.true_order())}, evaluations={result.evaluations}, "
        f"{result.seconds:.2f} s"
    )
    return 0


def _analyze(cfg: dict[str, Any], results: Path) -> int:
    """Statistics, tables and figures 4 to 11 from the raw runs in results/phase5."""
    import numpy as np

    from genomeasm.experiments import analysis, figures, runner

    summary = json.loads((results / "tuning_summary.json").read_text(encoding="utf-8"))
    runs, inst = analysis.load_runs(runner.read_records(results / "main_runs.jsonl.gz"))
    scaling_runs, scaling_inst = analysis.load_runs(
        runner.read_records(results / "scaling_runs.jsonl.gz")
    )
    best = summary["best_alns_policy"]
    tests = analysis.run_tests(runs, cfg, best)
    conv = analysis.convergence(
        runs, np.geomspace(0.02, cfg["budget_seconds"], 60), np.geomspace(1e3, 3e8, 60)
    )
    scaling = analysis.scaling_table(scaling_runs, scaling_inst, float(cfg["budget_seconds"]))
    by_algorithm = analysis.summarise(runs, ["algorithm"])
    by_instance = analysis.summarise(runs, ["instance", "algorithm"])
    true_ref = (
        inst.groupby("instance")[
            [c for c in inst.columns if c.startswith("bio_")] + ["n", "F_true"]
        ]
        .median()
        .reset_index()
    )
    out = results / "analysis"
    analysis.write_tables(
        {"by_algorithm": by_algorithm, "by_instance": by_instance, "tests": tests,
         "scaling": scaling, "true_order_reference": true_ref}, out)  # fmt: skip
    claim = analysis.claim_verdict(tests, float(cfg["statistics"]["alpha"]))
    info = {
        "best_alns_policy_from_tuning": best,
        "claim_bandit_beats_random": claim,
        "post_hoc_ucb1_tie_breaking": analysis.claim_verdict(
            tests, float(cfg["statistics"]["alpha"]), "post_hoc"
        ),
        "read_sets": int(runs.drop_duplicates(["instance", "seed_index"]).shape[0]),
        "provenance": provenance.collect(cfg, seed=int(cfg["seed"]), repo_dir=Path.cwd()),
    }
    _write_json(out / "analysis.json", info)
    _write_json(out / "convergence.json", conv)
    primary = tests[(tests["scope"] == "pooled") & (tests["family"] == "primary")]
    shown = {"by_algorithm": by_algorithm, "scaling": scaling, "tests_primary_pooled": primary}
    text = "\n\n".join(f"## {k}\n\n{analysis.markdown(v)}" for k, v in shown.items())
    (out / "tables.md").write_text(text, encoding="utf-8", newline="\n")
    paths = figures.make_figures(
        runs, tests, conv, scaling, float(cfg["budget_seconds"]), Path("reports/figures"),
        (results / "main_runs.jsonl.gz").as_posix(), info["provenance"]["config_hash"])  # fmt: skip
    print(f"wrote tables to {out} and {len(paths)} figure files; claim: {claim}")
    return 0


def _cmd_experiments(args: argparse.Namespace) -> int:
    """tune | run | scale: run a Phase 5 stage (resumable); analyze | figures are in analysis."""
    from genomeasm.experiments import plan, runner
    from genomeasm.experiments.design import load_config

    cfg = load_config(args.config)
    results = Path(args.results_dir)
    work = Path(args.work_dir)
    workers = args.workers or int(cfg["workers"])
    summary_path = results / "tuning_summary.json"
    if args.stage == "analyze":
        return _analyze(cfg, results)
    if args.stage == "tune":
        tasks = plan.tuning_tasks(cfg, seeds=args.seeds, trials=args.trials)
        runner.execute(tasks, work / "tuning.jsonl", workers)
        records = runner.read_records(work / "tuning.jsonl")
        frozen = (
            json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else None
        )
        summary = plan.select_tuned(cfg, records, frozen)
        stamp = provenance.collect(cfg, seed=int(cfg["seed"]), repo_dir=Path.cwd())
        if frozen:  # algorithms tuned earlier keep their entries; only the new ones are stamped
            summary["provenance"] = frozen["provenance"]
            added = dict(frozen.get("provenance_added", {}))
            added.update({a: stamp for a in summary["algorithms"] if a not in frozen["algorithms"]})
            summary["provenance_added"] = added
        else:
            summary["provenance"] = stamp
        _write_json(summary_path, summary)
        runner.compress(work / "tuning.jsonl", results / "tuning_runs.jsonl.gz")
        print(f"best ALNS policy on the tuning data: {summary['best_alns_policy']}")
        return 0
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if args.stage in ("rerun-lkh", "add-runs"):
        name = args.target
        records = runner.read_records(work / f"{name}.jsonl")
        algorithm = "lkh3" if args.stage == "rerun-lkh" else args.algorithm
        tasks = plan.missing_tasks(cfg, summary, records, algorithm, stage=name)
        runner.execute(tasks, work / f"{name}.jsonl", workers)
        runner.compress(work / f"{name}.jsonl", results / f"{name}_runs.jsonl.gz")
        return 0
    if args.stage == "run":
        specs = cfg["evaluation_instances"][: args.instances] if args.instances else None
        tasks = plan.main_tasks(cfg, summary, seeds=args.seeds, specs=specs)
        if args.seconds:
            for t in tasks:
                t["seconds"] = args.seconds
        name = "main"
    else:
        tasks = plan.scaling_tasks(cfg, summary, seeds=args.seeds)
        name = "scaling"
    runner.execute(tasks, work / f"{name}.jsonl", workers)
    runner.compress(work / f"{name}.jsonl", results / f"{name}_runs.jsonl.gz")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genomeasm", description=__doc__)
    parser.add_argument("--version", action="version", version=f"genomeasm {__version__}")
    top = parser.add_subparsers(dest="group", required=True)

    exp = top.add_parser("experiments", help="Phase 5 experiments (tune, run, scale)")
    exp.add_argument("stage", choices=["tune", "run", "scale", "analyze", "rerun-lkh", "add-runs"])
    exp.add_argument("--config", type=Path, default=Path("configs/phase5_experiments.yaml"))
    exp.add_argument("--results-dir", type=Path, default=Path("results/phase5"))
    exp.add_argument("--work-dir", type=Path, default=Path("data/processed/phase5_work"))
    exp.add_argument("--workers", type=int, default=0)
    exp.add_argument("--seeds", type=int, default=None, help="override the number of seeds")
    exp.add_argument("--trials", type=int, default=None, help="tune: override the trials")
    exp.add_argument("--instances", type=int, default=None, help="run: first k instances only")
    exp.add_argument("--seconds", type=float, default=None, help="run: override the budget")
    exp.add_argument(
        "--target", choices=["main", "scaling"], default="main", help="rerun-lkh, add-runs"
    )
    exp.add_argument("--algorithm", default="tabu_cl", help="add-runs: algorithm to add")
    exp.set_defaults(func=_cmd_experiments)

    solvers = top.add_parser("solvers", help="validate the solvers")
    solvers_sub = solvers.add_subparsers(dest="command", required=True)
    sv = solvers_sub.add_parser("validate", help="Phase 4 correctness validation")
    sv.add_argument("--config", type=Path, default=Path("configs/phase4_correctness.yaml"))
    sv.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    sv.add_argument("--results-dir", type=Path, default=Path("results/phase4"))
    sv.set_defaults(func=_cmd_solvers_validate)

    solve = top.add_parser("solve", help="run one solver on a simulated instance")
    solve.add_argument("accession")
    solve.add_argument("--coverage", type=float, required=True)
    solve.add_argument("--solver", choices=["greedy", "nn", "tabu", "alns", "lkh"], default="alns")
    solve.add_argument("--policy", default="ucb1", help="ALNS operator-selection policy")
    solve.add_argument("--seconds", type=float, default=5.0)
    solve.add_argument("--read-length", type=int, default=100)
    solve.add_argument("--sub-rate", type=float, default=0.0)
    solve.add_argument("--min-overlap", default="auto")
    solve.add_argument("--tolerance", type=float, default=0.0)
    solve.add_argument("--read-seed", type=int, default=0)
    solve.add_argument("--seed", type=int, default=0)
    solve.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    solve.set_defaults(func=_cmd_solve)

    docs = top.add_parser("docs", help="report, slides and oral-exam sheet (Phase 7)")
    docs_sub = docs.add_subparsers(dest="command", required=True)
    docs_build = docs_sub.add_parser("build", help="build all three from results/")
    docs_build.add_argument("--only", choices=["report", "slides", "oral"], default=None)
    docs_build.add_argument("--author", default="", help="name on the cover and the title slide")
    docs_build.set_defaults(func=_cmd_docs)
    docs_verify = docs_sub.add_parser(
        "verify", help="check the numbers in the documents and the charts"
    )
    docs_verify.add_argument("--out", type=Path, default=Path("results/phase8/verification.json"))
    docs_verify.set_defaults(func=_cmd_docs)

    app = top.add_parser("app", help="the interactive application (Phase 6)")
    app_sub = app.add_subparsers(dest="command", required=True)
    app_run = app_sub.add_parser("run", help="start the application (offline)")
    app_run.add_argument("--port", type=int, default=8501)
    app_run.add_argument("--browser", action="store_true", help="open it in the default browser")
    app_run.set_defaults(func=_cmd_app)
    app_bench = app_sub.add_parser("benchmark", help="measure the worst interactive runs")
    app_bench.add_argument(
        "--out", type=Path, default=Path("results/phase6/interactive_timing.json")
    )
    app_bench.set_defaults(func=_cmd_app)

    overlap = top.add_parser("overlap", help="validate the overlap graph against the ground truth")
    overlap_sub = overlap.add_subparsers(dest="command", required=True)
    ov_validate = overlap_sub.add_parser("validate", help="Phase 3 validation")
    ov_validate.add_argument("--config", type=Path, default=Path("configs/phase3_overlap.yaml"))
    ov_validate.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    ov_validate.add_argument("--results-dir", type=Path, default=Path("results/phase3"))
    ov_validate.add_argument("--figures-dir", type=Path, default=Path("reports/figures"))
    ov_validate.set_defaults(func=_cmd_overlap_validate)
    ov_figures = overlap_sub.add_parser("figures", help="redraw the Phase 3 figures")
    ov_figures.add_argument("--results-dir", type=Path, default=Path("results/phase3"))
    ov_figures.add_argument("--figures-dir", type=Path, default=Path("reports/figures"))
    ov_figures.set_defaults(func=_cmd_overlap_figures)

    asm = top.add_parser("assemble", help="walking skeleton: simulate, order, glue, vote")
    asm.add_argument("accession", help="versioned accession of a cached genome")
    asm.add_argument("--coverage", type=float, required=True)
    asm.add_argument("--read-length", type=int, default=100)
    asm.add_argument("--sub-rate", type=float, default=0.0)
    asm.add_argument("--min-overlap", default="auto", help="'auto' or a whole number")
    asm.add_argument("--tolerance", type=float, default=0.0, help="0 = exact overlaps")
    asm.add_argument("--order", choices=["trivial", "true", "random"], default="trivial")
    asm.add_argument("--seed", type=int, default=0)
    asm.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    asm.add_argument("--out-dir", type=Path, default=None)
    asm.set_defaults(func=_cmd_assemble)

    sim = top.add_parser("simulate", help="simulate sequencing reads and validate the simulator")
    sim_sub = sim.add_subparsers(dest="command", required=True)
    run = sim_sub.add_parser("run", help="simulate reads from a cached genome")
    run.add_argument("accession", help="versioned accession of a cached genome")
    run.add_argument("--coverage", type=float, required=True, help="nominal coverage c = N*l/G")
    run.add_argument("--read-length", type=int, default=100, help="mean read length")
    run.add_argument("--read-length-sd", type=float, default=0.0, help="0 = fixed length")
    run.add_argument("--sub-rate", type=float, default=0.0, help="substitution rate per base")
    run.add_argument("--ins-rate", type=float, default=0.0, help="insertion rate per base")
    run.add_argument("--del-rate", type=float, default=0.0, help="deletion rate per base")
    run.add_argument("--seed", type=int, default=0)
    run.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    run.add_argument("--out-dir", type=Path, default=None)
    run.set_defaults(func=_cmd_sim_run)

    validate = sim_sub.add_parser("validate", help="Phase 2 statistical validation")
    validate.add_argument("--config", type=Path, default=Path("configs/phase2_validation.yaml"))
    validate.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    validate.add_argument("--results-dir", type=Path, default=Path("results/phase2"))
    validate.add_argument("--figures-dir", type=Path, default=Path("reports/figures"))
    validate.set_defaults(func=_cmd_sim_validate)

    figures = sim_sub.add_parser("figures", help="redraw the Phase 2 figures from saved results")
    figures.add_argument("--results-dir", type=Path, default=Path("results/phase2"))
    figures.add_argument("--figures-dir", type=Path, default=Path("reports/figures"))
    figures.set_defaults(func=_cmd_sim_figures)

    data = top.add_parser("data", help="download, verify and describe genomes")
    sub = data.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("accession", nargs="*", help="versioned accessions (default: the config)")
    common.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    common.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)

    fetch = sub.add_parser("fetch", parents=[common], help="download genomes from NCBI")
    fetch.add_argument("--email", default=None, help="optional NCBI contact email (or NCBI_EMAIL)")
    fetch.add_argument("--force", action="store_true", help="download again even if cached")
    fetch.set_defaults(func=_cmd_fetch)

    verify = sub.add_parser("verify", parents=[common], help="check cached files (offline)")
    verify.set_defaults(func=_cmd_verify)

    cards = sub.add_parser("cards", parents=[common], help="write data cards (offline)")
    cards.add_argument("--processed-dir", type=Path, default=DEFAULT_PROCESSED_DIR)
    cards.set_defaults(func=_cmd_cards)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
