"""Genome object, alphabet rule (D-016), GC content, checksums and accession format."""

from __future__ import annotations

import hashlib

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genomeasm.io import (
    AccessionError,
    AmbiguousBaseError,
    FastaFormatError,
    Genome,
    format_fasta,
    parse_single_record,
    validate_accession,
)


def test_counts_length_and_gc() -> None:
    genome = Genome.from_sequence("NC_000001.1", "AACCGGTT")
    assert genome.length == 8
    assert genome.base_counts == {"A": 2, "C": 2, "G": 2, "T": 2}
    assert genome.gc_count == 4
    assert genome.gc_fraction == 0.5


def test_gc_extremes() -> None:
    assert Genome.from_sequence("NC_000001.1", "AATT").gc_fraction == 0.0
    assert Genome.from_sequence("NC_000001.1", "GGCC").gc_fraction == 1.0


@given(sequence=st.text(alphabet="ACGT", min_size=1, max_size=500))
def test_gc_matches_naive_count(sequence: str) -> None:
    genome = Genome.from_sequence("NC_000001.1", sequence)
    naive = sum(1 for c in sequence if c in "GC") / len(sequence)
    assert genome.gc_fraction == pytest.approx(naive)
    assert sum(genome.base_counts.values()) == genome.length == len(sequence)


def test_sequence_sha256_is_the_hash_of_the_letters() -> None:
    genome = Genome.from_sequence("NC_000001.1", "ACGTACGT")
    assert genome.sequence_sha256 == hashlib.sha256(b"ACGTACGT").hexdigest()


def test_sequence_sha256_ignores_line_wrapping_but_file_bytes_do_not() -> None:
    sequence = "ACGT" * 30
    narrow = format_fasta("NC_000001.1 x", sequence, width=10)
    wide = format_fasta("NC_000001.1 x", sequence, width=70)
    g_narrow = Genome.from_record(parse_single_record(narrow))
    g_wide = Genome.from_record(parse_single_record(wide))
    assert g_narrow.sequence_sha256 == g_wide.sequence_sha256
    assert hashlib.sha256(narrow.encode()).hexdigest() != hashlib.sha256(wide.encode()).hexdigest()


def test_ambiguity_code_is_rejected_with_position() -> None:
    with pytest.raises(AmbiguousBaseError, match=r"N at position 4") as info:
        Genome.from_sequence("NC_000001.1", "ACGTNACGT")
    assert "D-016" in str(info.value)


def test_several_ambiguity_codes_are_counted() -> None:
    with pytest.raises(AmbiguousBaseError, match="3 ambiguity symbol"):
        Genome.from_sequence("NC_000001.1", "ARCGYTN")


def test_ambiguity_in_a_parsed_record_is_rejected() -> None:
    record = parse_single_record(">NC_000001.1 x\nACGTNNNNACGT\n")
    assert "N" in record.sequence  # the parser keeps it ...
    with pytest.raises(AmbiguousBaseError):  # ... the genome refuses it
        Genome.from_record(record)


def test_other_invalid_symbols_are_a_format_error() -> None:
    with pytest.raises(FastaFormatError, match="invalid symbols"):
        Genome.from_sequence("NC_000001.1", "ACG1T")


def test_empty_sequence_is_rejected() -> None:
    with pytest.raises(FastaFormatError, match="empty"):
        Genome.from_sequence("NC_000001.1", "")


def test_from_record_checks_the_accession() -> None:
    record = parse_single_record(">NC_000001.1 x\nACGT\n")
    assert Genome.from_record(record, expected_accession="NC_000001.1").accession == "NC_000001.1"
    with pytest.raises(AccessionError, match="does not match"):
        Genome.from_record(record, expected_accession="NC_000002.1")


def test_lowercase_count_is_carried_into_the_genome() -> None:
    genome = Genome.from_record(parse_single_record(">NC_000001.1\nacgtAC\n"))
    assert genome.sequence == "ACGTAC"
    assert genome.n_lowercase_in_file == 4


@pytest.mark.parametrize("accession", ["NC_001422.1", "NC_001416.1", "AB123456.1", "U00096.3"])
def test_valid_accessions(accession: str) -> None:
    assert validate_accession(accession) == accession


@pytest.mark.parametrize(
    "accession",
    [
        "NC_001422",
        "nc_001422.1",
        "../x.1",
        "",
        "NC_001422.1; rm",
        "12345.1",
        "NC_001422.",
        "NC_1.1",
    ],
)
def test_invalid_accessions(accession: str) -> None:
    with pytest.raises(AccessionError):
        validate_accession(accession)
