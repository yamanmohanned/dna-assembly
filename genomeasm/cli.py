"""Command-line interface: ``genomeasm data fetch | verify | cards``."""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Sequence
from functools import partial
from pathlib import Path

import numpy as np
import yaml

from genomeasm import __version__, provenance
from genomeasm.datacard import write_cards
from genomeasm.io import (
    GenomeError,
    GenomeNotCachedError,
    download_genome,
    fetch_from_ncbi,
    load_genome,
)
from genomeasm.plots import make_phase2_figures
from genomeasm.sim_validation import run_validation
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads

DEFAULT_CONFIG = Path("configs/genomes.yaml")
DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_READS_DIR = Path("data/processed/reads")


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genomeasm", description=__doc__)
    parser.add_argument("--version", action="version", version=f"genomeasm {__version__}")
    top = parser.add_subparsers(dest="group", required=True)

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
