"""Download validation, the offline cache, checksums and the NCBI client (without network)."""

from __future__ import annotations

import hashlib
import io
import json
import urllib.error
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

from genomeasm.io import (
    AccessionError,
    AmbiguousBaseError,
    ChecksumMismatchError,
    GenomeNotCachedError,
    NcbiDataError,
    NcbiRecord,
    NetworkError,
    download_genome,
    fasta_path,
    fetch_from_ncbi,
    load_genome,
    meta_path,
    read_sums,
)

ACC = "NC_999999.1"
RecordFactory = Callable[..., NcbiRecord]


def _fixed_clock() -> datetime:
    return datetime(2026, 10, 22, 12, 0, 0, tzinfo=UTC)


class CountingFetcher:
    def __init__(self, record: NcbiRecord) -> None:
        self.record = record
        self.calls = 0

    def __call__(self, accession: str) -> NcbiRecord:
        self.calls += 1
        return self.record


def _no_fetch(accession: str) -> NcbiRecord:
    raise AssertionError("the fetcher must not be called")


# ---- download writes a complete, consistent cache -------------------------------------


def test_download_writes_fasta_meta_and_checksums(
    tmp_path: Path, ncbi_record: RecordFactory
) -> None:
    record = ncbi_record(ACC)
    loaded = download_genome(ACC, tmp_path, fetcher=lambda a: record, clock=_fixed_clock)

    fasta = fasta_path(tmp_path, ACC)
    assert fasta.read_bytes() == record.fasta_text.encode("utf-8")
    assert loaded.file_sha256 == hashlib.sha256(fasta.read_bytes()).hexdigest()
    assert read_sums(tmp_path) == {fasta.name: loaded.file_sha256}
    sums_line = (tmp_path / "SHA256SUMS").read_text().strip()
    assert sums_line == f"{loaded.file_sha256}  {fasta.name}"

    meta = json.loads(meta_path(tmp_path, ACC).read_text())
    assert meta["accession_version"] == ACC
    assert meta["ncbi_reported_length"] == loaded.genome.length
    assert meta["downloaded_utc"] == "2026-10-22T12:00:00+00:00"
    assert meta["esummary"]["accessionversion"] == ACC
    assert not list(tmp_path.glob("*.tmp"))


def test_two_genomes_keep_separate_sorted_checksum_lines(
    tmp_path: Path, ncbi_record: RecordFactory
) -> None:
    for accession in ("NC_000002.1", "NC_000001.1"):
        download_genome(accession, tmp_path, fetcher=CountingFetcher(ncbi_record(accession)))
    lines = (tmp_path / "SHA256SUMS").read_text().splitlines()
    assert [line.split("  ")[1] for line in lines] == ["NC_000001.1.fasta", "NC_000002.1.fasta"]


# ---- offline behaviour -----------------------------------------------------------------


def test_cached_genome_loads_with_the_network_blocked(
    tmp_path: Path, ncbi_record: RecordFactory, no_network: None
) -> None:
    download_genome(ACC, tmp_path, fetcher=lambda a: ncbi_record(ACC))
    loaded = load_genome(ACC, tmp_path, fetcher=_no_fetch)
    assert loaded.genome.accession == ACC
    assert loaded.genome.length == ncbi_record(ACC).summary["slen"]


def test_missing_genome_is_an_error_when_downloads_are_not_allowed(tmp_path: Path) -> None:
    with pytest.raises(GenomeNotCachedError, match="genomeasm data fetch"):
        load_genome(ACC, tmp_path, fetcher=_no_fetch)


def test_allowed_download_happens_once(tmp_path: Path, ncbi_record: RecordFactory) -> None:
    fetcher = CountingFetcher(ncbi_record(ACC))
    first = load_genome(ACC, tmp_path, allow_download=True, fetcher=fetcher)
    second = load_genome(ACC, tmp_path, allow_download=True, fetcher=fetcher)
    assert fetcher.calls == 1
    assert first.file_sha256 == second.file_sha256


def test_partial_cache_counts_as_not_cached(tmp_path: Path, ncbi_record: RecordFactory) -> None:
    download_genome(ACC, tmp_path, fetcher=lambda a: ncbi_record(ACC))
    meta_path(tmp_path, ACC).unlink()
    with pytest.raises(GenomeNotCachedError):
        load_genome(ACC, tmp_path)


# ---- tampering and inconsistent data ---------------------------------------------------


def test_modified_fasta_is_detected_by_the_checksum(
    tmp_path: Path, ncbi_record: RecordFactory
) -> None:
    download_genome(ACC, tmp_path, fetcher=lambda a: ncbi_record(ACC))
    path = fasta_path(tmp_path, ACC)
    text = path.read_text()
    first, rest = text.split("\n", 1)
    flipped = ("C" if rest[0] != "C" else "G") + rest[1:]
    path.write_text(first + "\n" + flipped)
    with pytest.raises(ChecksumMismatchError, match="SHA-256"):
        load_genome(ACC, tmp_path)


def test_malformed_checksum_file_is_rejected(tmp_path: Path, ncbi_record: RecordFactory) -> None:
    download_genome(ACC, tmp_path, fetcher=lambda a: ncbi_record(ACC))
    (tmp_path / "SHA256SUMS").write_text("not-a-digest  NC_999999.1.fasta\n")
    with pytest.raises(ChecksumMismatchError, match="malformed"):
        load_genome(ACC, tmp_path)


def test_meta_length_that_disagrees_with_the_fasta_is_rejected(
    tmp_path: Path, ncbi_record: RecordFactory
) -> None:
    download_genome(ACC, tmp_path, fetcher=lambda a: ncbi_record(ACC))
    path = meta_path(tmp_path, ACC)
    meta: dict[str, Any] = json.loads(path.read_text())
    meta["ncbi_reported_length"] += 1
    path.write_text(json.dumps(meta))
    with pytest.raises(NcbiDataError, match="cached NCBI length"):
        load_genome(ACC, tmp_path)


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"reported_length": 5}, NcbiDataError),
        ({"summary_accession": "NC_000002.1"}, NcbiDataError),
        ({"moltype": "rna"}, NcbiDataError),
        ({"header_accession": "NC_000002.1"}, AccessionError),
        ({"sequence": "ACGT" * 10 + "N" + "ACGT" * 10}, AmbiguousBaseError),
    ],
)
def test_bad_download_is_refused_and_nothing_is_cached(
    tmp_path: Path,
    ncbi_record: RecordFactory,
    kwargs: dict[str, Any],
    error: type[Exception],
) -> None:
    record = ncbi_record(ACC, **kwargs)
    with pytest.raises(error):
        download_genome(ACC, tmp_path, fetcher=lambda a: record)
    assert list(tmp_path.iterdir()) == []


def test_malformed_accession_never_touches_the_fetcher(tmp_path: Path) -> None:
    with pytest.raises(AccessionError):
        download_genome("../evil", tmp_path, fetcher=_no_fetch)


# ---- the NCBI client, with Entrez replaced by canned answers ---------------------------

_ESUMMARY = {
    "result": {
        "uids": ["42"],
        "42": {
            "uid": "42",
            "accessionversion": ACC,
            "organism": "Test phage",
            "moltype": "dna",
            "slen": 4,
        },
    }
}


def _patch_entrez(
    monkeypatch: pytest.MonkeyPatch, summary: object, fasta: object
) -> list[dict[str, Any]]:
    calls: list[dict[str, Any]] = []

    def esummary(**kwargs: Any) -> io.StringIO | io.BytesIO:
        calls.append({"fn": "esummary", **kwargs})
        if isinstance(summary, Exception):
            raise summary
        return io.BytesIO(json.dumps(summary).encode())  # NCBI sends JSON as bytes

    def efetch(**kwargs: Any) -> io.StringIO:
        calls.append({"fn": "efetch", **kwargs})
        if isinstance(fasta, Exception):
            raise fasta
        return io.StringIO(str(fasta))

    monkeypatch.setattr("Bio.Entrez.esummary", esummary)
    monkeypatch.setattr("Bio.Entrez.efetch", efetch)
    return calls


def test_fetch_from_ncbi_returns_record_and_provenance(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _patch_entrez(monkeypatch, _ESUMMARY, f">{ACC} Test phage\nACGT\n")
    record = fetch_from_ncbi(ACC)
    assert record.fasta_text == f">{ACC} Test phage\nACGT\n"
    assert record.summary["accessionversion"] == ACC
    assert record.provenance["email_sent"] == "no"
    assert record.provenance["entrez_tool"] == "genomeasm"
    assert [c["fn"] for c in calls] == ["esummary", "efetch"]
    assert calls[0]["db"] == "nuccore" and calls[1]["rettype"] == "fasta"


def test_fetch_from_ncbi_reports_when_an_email_was_sent(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_entrez(monkeypatch, _ESUMMARY, f">{ACC}\nACGT\n")
    assert fetch_from_ncbi(ACC, email="someone@example.org").provenance["email_sent"] == "yes"


def test_missing_email_warning_is_silenced_but_other_warnings_are_not(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import warnings

    biopython_text = (
        "\n            Email address is not specified.\n"
        "            To make use of NCBI's E-utilities, NCBI requires you to specify your\n"
    )

    def esummary(**kwargs: Any) -> io.BytesIO:
        warnings.warn(biopython_text, UserWarning, stacklevel=1)
        return io.BytesIO(json.dumps(_ESUMMARY).encode())

    def efetch(**kwargs: Any) -> io.StringIO:
        warnings.warn("some other problem", UserWarning, stacklevel=1)
        return io.StringIO(f">{ACC}\nACGT\n")

    monkeypatch.setattr("Bio.Entrez.esummary", esummary)
    monkeypatch.setattr("Bio.Entrez.efetch", efetch)
    with pytest.warns(UserWarning, match="some other problem") as caught:
        fetch_from_ncbi(ACC)
    assert not any("Email address" in str(w.message) for w in caught)


def test_unreachable_ncbi_is_a_network_error_and_invents_nothing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_entrez(monkeypatch, urllib.error.URLError("down"), "")
    with pytest.raises(NetworkError, match="No data is invented"):
        load_genome(ACC, tmp_path, allow_download=True)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "summary",
    [
        {"result": {"uids": []}},
        {"result": {"uids": ["1", "2"], "1": {}, "2": {}}},
        {"error": "bad"},
    ],
)
def test_unusable_ncbi_summary_is_rejected(
    monkeypatch: pytest.MonkeyPatch, summary: object
) -> None:
    _patch_entrez(monkeypatch, summary, f">{ACC}\nACGT\n")
    with pytest.raises(NcbiDataError):
        fetch_from_ncbi(ACC)
