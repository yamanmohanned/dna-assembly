"""Genome input: FASTA parsing, validation, NCBI download and an offline cache.

Plain-language story: a genome is one very long word written with four letters (A, C, G, T).
NCBI keeps the official copies. We download each copy once, check it carefully, save it on disk
together with a fingerprint (its SHA-256 checksum) and from then on only read the saved copy.
That way the program works without internet and always sees exactly the same letters.

Pipeline (full description with file names: docs/DATA_LOADING.md)::

    accession --Entrez--> FASTA text + NCBI summary --parse--> validate --write--> data/raw/
    data/raw/ --verify SHA-256--> parse --validate--> cross-check with NCBI summary --> Genome

Alphabet rule (decision D-016): a genome must consist only of the letters A, C, G, T.
Lowercase letters are accepted and turned into uppercase (the count is recorded). IUPAC
ambiguity codes (R Y S W K M B D H V N) are accepted by the FASTA parser, so that no
information is lost and they can be counted, but a ``Genome`` containing any of them is
rejected with ``AmbiguousBaseError``. Nothing is ever repaired, dropped or guessed.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import urllib.error
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

BASES = "ACGT"
AMBIGUITY_CODES = "RYSWKMBDHVN"
IUPAC_NUCLEOTIDE_CODES = frozenset(BASES + AMBIGUITY_CODES)
FASTA_WIDTH = 70
SUMS_FILENAME = "SHA256SUMS"

# Format check only: letters, optional underscore, 5 to 9 digits, a dot and a version number.
# Whether the accession really exists is checked against NCBI, not by this pattern.
_ACCESSION_RE = re.compile(r"[A-Z]{1,6}_?\d{5,9}\.\d+")
_SHA256_RE = re.compile(r"[0-9a-f]{64}")


class GenomeError(Exception):
    """Base class of every error raised by this module."""


class FastaFormatError(GenomeError, ValueError):
    """The text is not valid FASTA (empty, malformed, or contains invalid symbols)."""


class AmbiguousBaseError(GenomeError, ValueError):
    """The sequence contains symbols other than A, C, G, T (for example N)."""


class AccessionError(GenomeError, ValueError):
    """The accession is malformed or does not match the record that was read."""


class ChecksumMismatchError(GenomeError):
    """A cached file does not match its recorded SHA-256 checksum."""


class GenomeNotCachedError(GenomeError):
    """The genome is not in the local cache and downloading is not allowed."""


class NetworkError(GenomeError):
    """NCBI could not be reached. By rule R1 no data is ever invented in that case."""


class NcbiDataError(GenomeError):
    """NCBI (or the cache) returned data that contradicts itself or our expectations."""


# --------------------------------------------------------------------------------------
# FASTA parsing
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA entry: a header line and its sequence (uppercase, no whitespace)."""

    header: str
    sequence: str
    n_lowercase: int

    @property
    def identifier(self) -> str:
        """First word of the header (for NCBI records this is the versioned accession)."""
        return self.header.split(None, 1)[0]

    @property
    def description(self) -> str:
        """Everything after the identifier."""
        parts = self.header.split(None, 1)
        return parts[1] if len(parts) > 1 else ""


def _finish_record(header: str, header_line: int, chunks: list[str], n_lower: int) -> FastaRecord:
    if not chunks:
        raise FastaFormatError(f"line {header_line}: record {header!r} has no sequence")
    return FastaRecord(header=header, sequence="".join(chunks), n_lowercase=n_lower)


def parse_fasta(text: str) -> list[FastaRecord]:
    """Parse FASTA text into records.

    Accepted: Windows or Unix line endings, blank lines, lowercase letters, IUPAC ambiguity
    codes, several records. Rejected with ``FastaFormatError``: empty input, text before the
    first ``>``, an empty header, a header without sequence, non-ASCII characters, whitespace
    inside a sequence line, and any symbol outside the IUPAC nucleotide alphabet (digits,
    ``-``, ``*``, ``U`` and so on).
    """
    if text.startswith("﻿"):
        text = text[1:]
    records: list[FastaRecord] = []
    header: str | None = None
    header_line = 0
    chunks: list[str] = []
    n_lower = 0
    for line_no, raw in enumerate(text.replace("\r\n", "\n").replace("\r", "\n").split("\n"), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(_finish_record(header, header_line, chunks, n_lower))
            header = line[1:].strip()
            if not header:
                raise FastaFormatError(f"line {line_no}: empty header ('>' without identifier)")
            header_line, chunks, n_lower = line_no, [], 0
            continue
        if header is None:
            raise FastaFormatError(f"line {line_no}: sequence data before the first '>' header")
        if not line.isascii():
            raise FastaFormatError(f"line {line_no}: non-ASCII character in sequence")
        if any(ch.isspace() for ch in line):
            raise FastaFormatError(f"line {line_no}: whitespace inside a sequence line")
        upper = line.upper()
        bad = sorted(set(upper) - IUPAC_NUCLEOTIDE_CODES)
        if bad:
            raise FastaFormatError(f"line {line_no}: invalid nucleotide symbol(s) {bad}")
        n_lower += sum(ch.islower() for ch in line)
        chunks.append(upper)
    if header is not None:
        records.append(_finish_record(header, header_line, chunks, n_lower))
    if not records:
        raise FastaFormatError("no FASTA record found (empty input)")
    return records


def parse_single_record(text: str) -> FastaRecord:
    """Parse FASTA text that must contain exactly one record (one genome per file)."""
    records = parse_fasta(text)
    if len(records) != 1:
        raise FastaFormatError(f"expected exactly one FASTA record, found {len(records)}")
    return records[0]


def format_fasta(header: str, sequence: str, width: int = FASTA_WIDTH) -> str:
    """Write one record as FASTA text with ``width`` letters per line."""
    if width < 1:
        raise ValueError("width must be at least 1")
    lines = [f">{header}"]
    lines += [sequence[i : i + width] for i in range(0, len(sequence), width)]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------------------
# Genome
# --------------------------------------------------------------------------------------


def validate_accession(accession: str) -> str:
    """Return ``accession`` if it looks like a versioned accession such as ``NC_001422.1``."""
    if not _ACCESSION_RE.fullmatch(accession):
        raise AccessionError(
            f"{accession!r} is not a versioned accession (expected something like NC_001422.1)"
        )
    return accession


@dataclass(frozen=True)
class Genome:
    """A validated genome: uppercase A/C/G/T only, with fingerprints and base counts.

    Real-data rule R1: genomes used for results come from NCBI. ``from_sequence`` also lets
    unit tests build tiny synthetic strings; they must never be used for reported results.
    """

    accession: str
    description: str
    sequence: str
    base_counts: dict[str, int]
    n_lowercase_in_file: int = 0

    @classmethod
    def from_sequence(
        cls,
        accession: str,
        sequence: str,
        *,
        description: str = "",
        n_lowercase_in_file: int = 0,
    ) -> Genome:
        """Validate an uppercase sequence and build a ``Genome``."""
        if not sequence:
            raise FastaFormatError(f"{accession}: empty sequence")
        counts = {base: sequence.count(base) for base in BASES}
        n_other = len(sequence) - sum(counts.values())
        if n_other:
            offenders = [(i, c) for i, c in enumerate(sequence) if c not in BASES]
            symbols = sorted({c for _, c in offenders})
            first = ", ".join(f"{c} at position {i}" for i, c in offenders[:5])
            if set(symbols) <= set(AMBIGUITY_CODES):
                raise AmbiguousBaseError(
                    f"{accession}: {n_other} ambiguity symbol(s) {symbols} (first: {first}). "
                    "Rule D-016: genomes must be pure A/C/G/T; nothing is repaired or dropped."
                )
            raise FastaFormatError(f"{accession}: invalid symbols {symbols} (first: {first})")
        return cls(
            accession=accession,
            description=description,
            sequence=sequence,
            base_counts=counts,
            n_lowercase_in_file=n_lowercase_in_file,
        )

    @classmethod
    def from_record(cls, record: FastaRecord, *, expected_accession: str | None = None) -> Genome:
        """Build a ``Genome`` from a parsed FASTA record, checking the accession if given."""
        if expected_accession is not None and record.identifier != expected_accession:
            raise AccessionError(
                f"FASTA identifier {record.identifier!r} does not match the requested "
                f"accession {expected_accession!r}"
            )
        return cls.from_sequence(
            record.identifier,
            record.sequence,
            description=record.description,
            n_lowercase_in_file=record.n_lowercase,
        )

    @property
    def length(self) -> int:
        """Number of letters (base pairs for a double-stranded genome)."""
        return len(self.sequence)

    @property
    def gc_count(self) -> int:
        """Number of G and C letters."""
        return self.base_counts["G"] + self.base_counts["C"]

    @property
    def gc_fraction(self) -> float:
        """GC content: (G + C) / length, a number between 0 and 1."""
        return self.gc_count / self.length

    @property
    def sequence_sha256(self) -> str:
        """SHA-256 of the uppercase sequence letters only (independent of line wrapping)."""
        return hashlib.sha256(self.sequence.encode("ascii")).hexdigest()


# --------------------------------------------------------------------------------------
# NCBI download
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class NcbiRecord:
    """What NCBI returned for one accession: FASTA text plus its own summary of the record."""

    fasta_text: str
    summary: dict[str, Any]
    provenance: dict[str, str]


Fetcher = Callable[[str], NcbiRecord]


def _as_text(data: str | bytes) -> str:
    return data.decode("utf-8") if isinstance(data, bytes) else data


def fetch_from_ncbi(
    accession: str,
    *,
    tool: str = "genomeasm",
    email: str | None = None,
    api_key: str | None = None,
) -> NcbiRecord:
    """Download one record from NCBI nuccore with Entrez (``esummary`` then ``efetch``).

    Biopython limits the request rate to 3 per second without an API key. The contact email
    is optional; without one Biopython would print a warning, which is silenced here on
    purpose because sending no email is a deliberate choice (decision D-010).
    """
    import Bio
    from Bio import Entrez

    validate_accession(accession)
    Entrez.tool = tool
    Entrez.email = email  # type: ignore[assignment]  # Biopython declares these as None
    Entrez.api_key = api_key  # type: ignore[assignment]
    try:
        with warnings.catch_warnings():
            if email is None:
                # (?s): Biopython's message starts with a newline, which "." would not match.
                warnings.filterwarnings(
                    "ignore", message=r"(?s).*Email address is not specified", category=UserWarning
                )
            with Entrez.esummary(db="nuccore", id=accession, retmode="json") as handle:
                payload = _as_text(handle.read())
            with Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text") as h:
                fasta_text = _as_text(h.read())
    except (urllib.error.URLError, OSError, TimeoutError) as exc:
        raise NetworkError(
            f"cannot reach NCBI for {accession}: {exc}. No data is invented (rule R1)."
        ) from exc
    try:
        result = json.loads(payload)["result"]
        uids = result["uids"]
        if len(uids) != 1:
            raise NcbiDataError(f"NCBI returned {len(uids)} records for {accession}")
        summary = dict(result[uids[0]])
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        raise NcbiDataError(f"unreadable NCBI summary for {accession}: {exc!r}") from exc
    return NcbiRecord(
        fasta_text=fasta_text,
        summary=summary,
        provenance={
            "entrez_tool": tool,
            "biopython_version": Bio.__version__,
            "database": "nuccore",
            "email_sent": "yes" if email else "no",
        },
    )


# --------------------------------------------------------------------------------------
# Offline cache
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class LoadedGenome:
    """A genome read from the cache, with the NCBI metadata and the file fingerprint."""

    genome: Genome
    meta: dict[str, Any]
    file_sha256: str
    fasta_path: Path


def fasta_path(raw_dir: Path, accession: str) -> Path:
    return Path(raw_dir) / f"{validate_accession(accession)}.fasta"


def meta_path(raw_dir: Path, accession: str) -> Path:
    return Path(raw_dir) / f"{validate_accession(accession)}.meta.json"


def sha256_file(path: Path) -> str:
    """SHA-256 of the exact bytes of a file."""
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_sums(raw_dir: Path) -> dict[str, str]:
    """Read ``SHA256SUMS`` (the format understood by ``sha256sum -c``) into name -> digest."""
    path = Path(raw_dir) / SUMS_FILENAME
    sums: dict[str, str] = {}
    if not path.is_file():
        return sums
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        digest, _, name = line.partition(" ")
        name = name.lstrip(" *")
        if not _SHA256_RE.fullmatch(digest) or not name:
            raise ChecksumMismatchError(f"{path.name} line {number} is malformed")
        sums[name] = digest
    return sums


def _atomic_write(path: Path, data: bytes) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _write_sums(raw_dir: Path, sums: dict[str, str]) -> None:
    text = "".join(f"{digest}  {name}\n" for name, digest in sorted(sums.items()))
    _atomic_write(Path(raw_dir) / SUMS_FILENAME, text.encode("utf-8"))


def download_genome(
    accession: str,
    raw_dir: Path,
    *,
    fetcher: Fetcher | None = None,
    clock: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> LoadedGenome:
    """Download, validate and cache one genome. Nothing is written unless every check passes.

    Checks: FASTA parses, identifier equals the requested accession, sequence is pure A/C/G/T,
    NCBI's own summary agrees on accession, molecule type (dna) and length.
    """
    validate_accession(accession)
    record = (fetcher or fetch_from_ncbi)(accession)
    genome = Genome.from_record(
        parse_single_record(record.fasta_text), expected_accession=accession
    )
    summary = record.summary
    if summary.get("accessionversion") != accession:
        raise NcbiDataError(
            f"NCBI summary is for {summary.get('accessionversion')!r}, not {accession!r}"
        )
    if summary.get("moltype") != "dna":
        raise NcbiDataError(f"{accession}: molecule type {summary.get('moltype')!r}, expected dna")
    try:
        reported_length = int(summary["slen"])
    except (KeyError, TypeError, ValueError) as exc:
        raise NcbiDataError(f"{accession}: NCBI summary has no usable length") from exc
    if reported_length != genome.length:
        raise NcbiDataError(
            f"{accession}: NCBI reports {reported_length} letters, FASTA has {genome.length}"
        )

    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    data = record.fasta_text.encode("utf-8")
    meta = {
        "accession_version": accession,
        "ncbi_uid": summary.get("uid"),
        "organism": summary.get("organism"),
        "title": summary.get("title"),
        "moltype": summary.get("moltype"),
        "ncbi_reported_length": reported_length,
        "source_database": summary.get("sourcedb"),
        "completeness": summary.get("completeness"),
        "source": "NCBI nuccore via Entrez (esummary + efetch)",
        "downloaded_utc": clock().isoformat(timespec="seconds"),
        "provenance": dict(record.provenance),
        "esummary": summary,
    }
    fasta_file = fasta_path(raw, accession)
    _atomic_write(fasta_file, data)
    _atomic_write(
        meta_path(raw, accession),
        (json.dumps(meta, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    sums = read_sums(raw)
    sums[fasta_file.name] = hashlib.sha256(data).hexdigest()
    _write_sums(raw, sums)
    return load_genome(accession, raw)


def load_genome(
    accession: str,
    raw_dir: Path,
    *,
    allow_download: bool = False,
    fetcher: Fetcher | None = None,
) -> LoadedGenome:
    """Load a genome from the offline cache, verifying its checksum first.

    By default nothing is downloaded (the application must work offline). Set
    ``allow_download=True`` to fetch a missing genome from NCBI.
    """
    validate_accession(accession)
    raw = Path(raw_dir)
    fasta_file, meta_file = fasta_path(raw, accession), meta_path(raw, accession)
    sums = read_sums(raw)
    expected_digest = sums.get(fasta_file.name)
    if not (fasta_file.is_file() and meta_file.is_file() and expected_digest):
        if not allow_download:
            raise GenomeNotCachedError(
                f"{accession} is not fully cached in {raw} (needs {fasta_file.name}, "
                f"{meta_file.name} and a {SUMS_FILENAME} entry). Run 'genomeasm data fetch'."
            )
        return download_genome(accession, raw, fetcher=fetcher)

    actual_digest = sha256_file(fasta_file)
    if actual_digest != expected_digest:
        raise ChecksumMismatchError(
            f"{fasta_file.name}: SHA-256 is {actual_digest}, {SUMS_FILENAME} says {expected_digest}"
        )
    genome = Genome.from_record(
        parse_single_record(fasta_file.read_text(encoding="utf-8")),
        expected_accession=accession,
    )
    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    if meta.get("accession_version") != accession:
        raise NcbiDataError(f"{meta_file.name} belongs to {meta.get('accession_version')!r}")
    if meta.get("ncbi_reported_length") != genome.length:
        raise NcbiDataError(
            f"{accession}: cached NCBI length {meta.get('ncbi_reported_length')} "
            f"!= parsed length {genome.length}"
        )
    return LoadedGenome(genome=genome, meta=meta, file_sha256=actual_digest, fasta_path=fasta_file)
