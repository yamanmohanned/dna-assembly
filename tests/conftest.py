"""Shared test helpers. Synthetic sequences here are for unit tests only, never for results."""

from __future__ import annotations

import random
import socket
from collections.abc import Callable
from typing import Any

import pytest

from genomeasm.io import NcbiRecord, format_fasta

_rng = random.Random(20261022)
DEFAULT_SEQUENCE = "".join(_rng.choice("ACGT") for _ in range(200))


@pytest.fixture
def ncbi_record() -> Callable[..., NcbiRecord]:
    """Factory for fake NCBI answers; every field can be corrupted to test the checks."""

    def factory(
        accession: str = "NC_999999.1",
        sequence: str | None = None,
        *,
        reported_length: int | None = None,
        moltype: str = "dna",
        summary_accession: str | None = None,
        header_accession: str | None = None,
        organism: str = "Test phage",
    ) -> NcbiRecord:
        seq = DEFAULT_SEQUENCE if sequence is None else sequence
        header = f"{header_accession or accession} {organism}, complete genome"
        summary: dict[str, Any] = {
            "uid": "1",
            "accessionversion": summary_accession or accession,
            "organism": organism,
            "title": f"{organism}, complete genome",
            "moltype": moltype,
            "slen": len(seq) if reported_length is None else reported_length,
            "sourcedb": "refseq",
            "completeness": "complete",
            "strand": "ss",
            "updatedate": "2023/01/11",
            "topology": "circular",
        }
        return NcbiRecord(
            fasta_text=format_fasta(header, seq),
            summary=summary,
            provenance={"entrez_tool": "test"},
        )

    return factory


@pytest.fixture
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make any attempt to use the network fail loudly."""

    def boom(*args: object, **kwargs: object) -> None:
        raise AssertionError("network access attempted while the test requires offline mode")

    monkeypatch.setattr(socket, "socket", boom)
    monkeypatch.setattr(socket, "create_connection", boom)
    monkeypatch.setattr(socket, "getaddrinfo", boom)
