"""FASTA parser: valid, empty, malformed, lowercase and ambiguity-code inputs."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genomeasm.io import (
    FastaFormatError,
    format_fasta,
    parse_fasta,
    parse_single_record,
)


def test_valid_single_record() -> None:
    record = parse_single_record(">NC_001422.1 Test phage, complete genome\nACGT\nTTGA\n")
    assert record.header == "NC_001422.1 Test phage, complete genome"
    assert record.identifier == "NC_001422.1"
    assert record.description == "Test phage, complete genome"
    assert record.sequence == "ACGTTTGA"
    assert record.n_lowercase == 0


def test_header_without_description() -> None:
    record = parse_single_record(">only_id\nACGT\n")
    assert record.identifier == "only_id"
    assert record.description == ""


@pytest.mark.parametrize("newline", ["\n", "\r\n", "\r"])
def test_line_endings_are_equivalent(newline: str) -> None:
    text = newline.join([">id.1 d", "ACGT", "TTGA", ""])
    assert parse_single_record(text).sequence == "ACGTTTGA"


def test_blank_lines_and_trailing_spaces_are_ignored() -> None:
    record = parse_single_record("\n>id.1 d  \n\nACGT  \n\n  TTGA\n\n\n")
    assert record.sequence == "ACGTTTGA"


def test_byte_order_mark_is_ignored() -> None:
    assert parse_single_record("﻿>id.1\nACGT\n").sequence == "ACGT"


def test_lowercase_is_uppercased_and_counted() -> None:
    record = parse_single_record(">id.1\nacgtAC\nggTT\n")
    assert record.sequence == "ACGTACGGTT"
    assert record.n_lowercase == 6


def test_ambiguity_codes_are_kept_by_the_parser() -> None:
    record = parse_single_record(">id.1\nACGTNRYSWKMBDHVn\n")
    assert record.sequence == "ACGTNRYSWKMBDHVN"
    assert record.n_lowercase == 1


def test_several_records() -> None:
    records = parse_fasta(">a.1\nAC\n>b.1\nGT\n")
    assert [(r.identifier, r.sequence) for r in records] == [("a.1", "AC"), ("b.1", "GT")]


def test_single_record_parser_rejects_two_records() -> None:
    with pytest.raises(FastaFormatError, match="exactly one"):
        parse_single_record(">a.1\nAC\n>b.1\nGT\n")


@pytest.mark.parametrize("text", ["", "   \n\n\t\n", "﻿"])
def test_empty_input_is_rejected(text: str) -> None:
    with pytest.raises(FastaFormatError, match="empty input"):
        parse_fasta(text)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (">id.1\n", "no sequence"),
        (">id.1\n>id.2\nACGT\n", "no sequence"),
        (">\nACGT\n", "empty header"),
        ("ACGT\n>id.1\nACGT\n", "before the first"),
        (">id.1\nAC1T\n", "invalid nucleotide"),
        (">id.1\nACG-T\n", "invalid nucleotide"),
        (">id.1\nACGU\n", "invalid nucleotide"),
        (">id.1\nAC*T\n", "invalid nucleotide"),
        (">id.1\nACG T\n", "whitespace"),
        (">id.1\nAC\tGT\n", "whitespace"),
        (">id.1\nACGß\n", "non-ASCII"),
        (">id.1\n;comment\nACGT\n", "invalid nucleotide"),
    ],
)
def test_malformed_input_is_rejected(text: str, message: str) -> None:
    with pytest.raises(FastaFormatError, match=message):
        parse_fasta(text)


def test_error_message_names_the_line() -> None:
    with pytest.raises(FastaFormatError, match="line 3"):
        parse_fasta(">id.1\nACGT\nAC9T\n")


def test_format_fasta_wraps_lines() -> None:
    assert format_fasta("id.1 d", "ACGTACGT", width=3) == ">id.1 d\nACG\nTAC\nGT\n"


def test_format_fasta_rejects_bad_width() -> None:
    with pytest.raises(ValueError, match="width"):
        format_fasta("id.1", "ACGT", width=0)


@given(
    sequence=st.text(alphabet="ACGT", min_size=1, max_size=400),
    width=st.integers(min_value=1, max_value=120),
)
def test_format_then_parse_roundtrip(sequence: str, width: int) -> None:
    record = parse_single_record(format_fasta("id.1 some description", sequence, width=width))
    assert record.sequence == sequence
    assert record.header == "id.1 some description"
    assert record.n_lowercase == 0


@given(sequence=st.text(alphabet="ACGTacgt", min_size=1, max_size=300))
def test_case_folding_property(sequence: str) -> None:
    record = parse_single_record(format_fasta("id.1", sequence))
    assert record.sequence == sequence.upper()
    assert record.n_lowercase == sum(c.islower() for c in sequence)
