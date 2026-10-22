"""Command-line interface: ``genomeasm data fetch | verify | cards``."""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence
from functools import partial
from pathlib import Path

import yaml

from genomeasm import __version__
from genomeasm.datacard import write_cards
from genomeasm.io import (
    GenomeError,
    GenomeNotCachedError,
    download_genome,
    fetch_from_ncbi,
    load_genome,
)

DEFAULT_CONFIG = Path("configs/genomes.yaml")
DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PROCESSED_DIR = Path("data/processed")


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genomeasm", description=__doc__)
    parser.add_argument("--version", action="version", version=f"genomeasm {__version__}")
    top = parser.add_subparsers(dest="group", required=True)
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
