"""Checks on the real NCBI files committed in data/raw and the data cards in data/processed.

These tests run with the network blocked: they prove that the saved copies alone are enough.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

from genomeasm import cli
from genomeasm.datacard import build_card, render_card_markdown
from genomeasm.io import load_genome, read_sums

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
ACCESSIONS = cli.load_accessions(ROOT / "configs" / "genomes.yaml")

pytestmark = pytest.mark.usefixtures("no_network")


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_checksum_file_matches_the_bytes_on_disk(accession: str) -> None:
    data = (RAW / f"{accession}.fasta").read_bytes()
    assert read_sums(RAW)[f"{accession}.fasta"] == hashlib.sha256(data).hexdigest()


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_genome_loads_offline_and_agrees_with_ncbi(accession: str) -> None:
    loaded = load_genome(accession, RAW)
    assert loaded.genome.accession == accession
    assert loaded.genome.length == loaded.meta["ncbi_reported_length"]
    assert loaded.meta["moltype"] == "dna"
    assert loaded.meta["source_database"] == "refseq"
    assert loaded.genome.n_lowercase_in_file == 0


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_independent_recount_of_the_raw_file(accession: str) -> None:
    """Recount the letters without any genomeasm code and compare with the Genome object."""
    lines = (RAW / f"{accession}.fasta").read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith(f">{accession} ")
    letters = Counter("".join(line.strip() for line in lines[1:]))
    assert set(letters) == set("ACGT"), f"unexpected symbols: {sorted(letters)}"

    genome = load_genome(accession, RAW).genome
    assert dict(letters) == genome.base_counts
    assert sum(letters.values()) == genome.length
    assert (letters["G"] + letters["C"]) / sum(letters.values()) == genome.gc_fraction


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_data_cards_on_disk_were_generated_from_the_cached_files(accession: str) -> None:
    card = build_card(load_genome(accession, RAW))
    stored = json.loads((PROCESSED / f"{accession}.card.json").read_text(encoding="utf-8"))
    assert stored == card
    assert (PROCESSED / f"{accession}.card.md").read_text(encoding="utf-8") == render_card_markdown(
        card
    )


def test_every_data_card_belongs_to_a_configured_genome() -> None:
    on_disk = {p.name.removesuffix(".card.json") for p in PROCESSED.glob("*.card.json")}
    assert on_disk == set(ACCESSIONS)
