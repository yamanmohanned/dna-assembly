"""Data cards, the command-line interface and the glossary file."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pytest
import yaml

from genomeasm import cli
from genomeasm.datacard import build_card, render_card_markdown, write_cards
from genomeasm.io import (
    NcbiRecord,
    NetworkError,
    download_genome,
    load_genome,
    validate_accession,
)

ACC = "NC_999999.1"
ROOT = Path(__file__).resolve().parent.parent
RecordFactory = Callable[..., NcbiRecord]


def _cached(tmp_path: Path, ncbi_record: RecordFactory, accession: str = ACC) -> Path:
    record = ncbi_record(accession, sequence="AACCGGTT" * 5)
    download_genome(accession, tmp_path, fetcher=lambda a: record)
    return tmp_path


def test_card_contains_the_required_facts(tmp_path: Path, ncbi_record: RecordFactory) -> None:
    card = build_card(load_genome(ACC, _cached(tmp_path, ncbi_record)))
    assert card["accession_version"] == ACC
    assert card["accession"] == "NC_999999"
    assert card["version"] == 1
    assert card["organism"] == "Test phage"
    assert card["length_bp"] == 40
    assert card["length_matches_ncbi"] is True
    assert card["gc_fraction"] == 0.5
    assert card["gc_percent"] == 50.0
    assert card["base_counts"] == {"A": 10, "C": 10, "G": 10, "T": 10}
    assert card["non_acgt_symbols"] == 0
    assert card["ncbi_source_database"] == "refseq"
    assert card["ncbi_completeness"] == "complete"
    assert card["topology_reported_by_ncbi"] == "circular"
    assert card["strand_reported_by_ncbi"] == "ss"
    assert len(card["file_sha256"]) == 64 and len(card["sequence_sha256"]) == 64


def test_cards_are_deterministic_and_use_only_the_cache(
    tmp_path: Path, ncbi_record: RecordFactory, no_network: None
) -> None:
    raw = _cached(tmp_path / "raw", ncbi_record)
    first = write_cards(raw, tmp_path / "out", [ACC])
    snapshot = {p.name: p.read_bytes() for p in first}
    second = write_cards(raw, tmp_path / "out", [ACC])
    assert {p.name: p.read_bytes() for p in second} == snapshot
    assert sorted(snapshot) == [f"{ACC}.card.json", f"{ACC}.card.md"]


def test_markdown_card_shows_gc_and_checksum(tmp_path: Path, ncbi_record: RecordFactory) -> None:
    card = build_card(load_genome(ACC, _cached(tmp_path, ncbi_record)))
    text = render_card_markdown(card)
    assert f"# Data card: {ACC}" in text
    assert "50.00 %" in text
    assert card["file_sha256"] in text
    assert "40 bp" in text


# ---- command line ------------------------------------------------------------------------


def test_cli_verify_cards_and_fetch_on_a_cached_genome(
    tmp_path: Path,
    ncbi_record: RecordFactory,
    capsys: pytest.CaptureFixture[str],
    no_network: None,
) -> None:
    raw = _cached(tmp_path / "raw", ncbi_record)
    base = ["data"]
    assert cli.main([*base, "verify", ACC, "--raw-dir", str(raw)]) == 0
    assert f"{ACC}: OK, 40 bp" in capsys.readouterr().out

    out = tmp_path / "processed"
    assert cli.main([*base, "cards", ACC, "--raw-dir", str(raw), "--processed-dir", str(out)]) == 0
    assert (out / f"{ACC}.card.json").is_file()

    assert cli.main([*base, "fetch", ACC, "--raw-dir", str(raw)]) == 0
    assert "already cached and verified" in capsys.readouterr().out


def test_cli_verify_reports_a_missing_genome(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert cli.main(["data", "verify", ACC, "--raw-dir", str(tmp_path)]) == 1
    assert "error:" in capsys.readouterr().err


def test_cli_fetch_reports_network_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def down(accession: str, **kwargs: object) -> NcbiRecord:
        raise NetworkError("cannot reach NCBI")

    monkeypatch.setattr(cli, "fetch_from_ncbi", down)
    assert cli.main(["data", "fetch", ACC, "--raw-dir", str(tmp_path)]) == 1
    assert "cannot reach NCBI" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_cli_help_exits_cleanly() -> None:
    with pytest.raises(SystemExit) as info:
        cli.main(["data", "fetch", "--help"])
    assert info.value.code == 0


def test_module_entry_point_runs() -> None:
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "genomeasm", "--version"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=False,
    )
    assert result.returncode == 0
    assert result.stdout.startswith("genomeasm ")


def test_project_config_lists_valid_accessions() -> None:
    accessions = cli.load_accessions(ROOT / "configs" / "genomes.yaml")
    assert len(accessions) >= 2
    assert [validate_accession(a) for a in accessions] == accessions
    assert len(set(accessions)) == len(accessions)


# ---- glossary ----------------------------------------------------------------------------


def test_glossary_entries_are_complete_and_unique() -> None:
    data = yaml.safe_load((ROOT / "docs" / "glossary.yaml").read_text(encoding="utf-8"))
    terms = data["terms"]
    assert terms, "glossary is empty"
    for entry in terms:
        for key in ("term", "plain", "precise", "phase"):
            assert entry.get(key), f"{entry.get('term')!r} lacks {key}"
    names = [entry["term"].lower() for entry in terms]
    assert len(names) == len(set(names)), "duplicate glossary terms"


def test_json_card_roundtrips(tmp_path: Path, ncbi_record: RecordFactory) -> None:
    raw = _cached(tmp_path / "raw", ncbi_record)
    paths = write_cards(raw, tmp_path / "out", [ACC])
    stored = json.loads(paths[0].read_text())
    assert stored == build_card(load_genome(ACC, raw))
