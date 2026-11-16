"""Knowledge checks of the molecular-biology pages.

Wherever a question has a computable answer, the correct option is computed by the same functions
the application uses (complement, transcription, the codon table, the mutation classifier), so a
question cannot disagree with the rest of the application. A test checks every question.
"""

from __future__ import annotations

from dataclasses import dataclass

from genomeasm.molbio import sequence

EXAMPLE = "ATGTCTAAAGGT"  # the first 12 letters of gene J of phiX174


@dataclass(frozen=True)
class Question:
    topic: str
    text: str
    options: tuple[str, ...]
    answer: int  # index of the correct option
    why: str


def _place(correct: str, others: tuple[str, ...], position: int) -> tuple[tuple[str, ...], int]:
    """Options with the correct one at a fixed position (so the order is the same on every run)."""
    options = list(others)
    options.insert(position, correct)
    if len(set(options)) != len(options):
        raise ValueError("two options of a question are identical")
    return tuple(options), position


def questions() -> list[Question]:
    q = []
    partner = sequence.reverse_complement(EXAMPLE)
    options, answer = _place(
        partner, (EXAMPLE[::-1], sequence.complement(EXAMPLE), sequence.transcribe(EXAMPLE)), 2
    )
    q.append(Question(
        "structure",
        f"One strand reads 5′-{EXAMPLE}-3′. What is the other strand, also written 5′ → 3′?",
        tuple(f"5′-{o}-3′" for o in options), answer,
        "Each base is replaced by its partner (A↔T, G↔C) and, because the strands are "
        "antiparallel, the result is read in the opposite direction: the reverse complement.",
    ))  # fmt: skip
    bonds = sequence.hydrogen_bonds(EXAMPLE)
    options, answer = _place(
        str(bonds), (str(2 * len(EXAMPLE)), str(3 * len(EXAMPLE)), str(len(EXAMPLE))), 1
    )
    n_gc = EXAMPLE.count("G") + EXAMPLE.count("C")
    q.append(Question(
        "structure",
        f"How many Watson-Crick hydrogen bonds hold {EXAMPLE} to its partner strand?",
        options, answer,
        f"2 for each A·T pair and 3 for each G·C pair: {len(EXAMPLE) - n_gc} × 2 + {n_gc} × 3 "
        f"= {bonds}.",
    ))  # fmt: skip
    q.append(Question(
        "structure",
        "Which parts of a nucleotide form the backbone of a DNA strand?",
        ("The bases", "Sugar and phosphate", "Phosphate and base", "Hydrogen bonds"), 1,
        "Sugars and phosphates alternate along the strand, joined by phosphodiester bonds from "
        "the 3′ carbon of one sugar to the 5′ carbon of the next. The bases point inwards.",
    ))  # fmt: skip
    q.append(Question(
        "structure",
        "A genome contains 24% A and 31% T. What does this suggest?",
        ("It is double-stranded DNA", "It is single-stranded", "It contains no G or C",
         "The measurement must be wrong"), 1,
        "In a double helix every A faces a T, so their amounts are equal. Unequal amounts mean "
        "there is no partner strand. This is the case of phiX174 inside the virus particle.",
    ))  # fmt: skip
    q.append(Question(
        "replication",
        "Why is one new strand made in short Okazaki fragments?",
        ("DNA polymerase is too slow", "DNA polymerase only extends 3′ ends, and the two "
         "templates are antiparallel", "The helicase cuts the template",
         "Primers can only be made on one strand"), 1,
        "Synthesis always runs 5′ → 3′. On one template this follows the fork; on the other it "
        "runs away from the fork, so it must restart again and again.",
    ))  # fmt: skip
    q.append(Question(
        "replication",
        "After one round of replication, each daughter double helix contains …",
        ("two new strands", "two old strands", "one old and one new strand",
         "pieces of old and new mixed in both strands"), 2,
        "Replication is semiconservative: each old strand is the template for one new strand "
        "and stays paired with it (Meselson and Stahl 1958).",
    ))  # fmt: skip
    rna = sequence.transcribe(EXAMPLE)
    template = sequence.complement(EXAMPLE)
    options, answer = _place(
        rna, (sequence.transcribe(template), EXAMPLE, sequence.transcribe(EXAMPLE[::-1])), 0
    )
    q.append(Question(
        "expression",
        f"The coding strand of a gene begins 5′-{EXAMPLE}-3′. How does its mRNA begin?",
        tuple(f"5′-{o}-3′" for o in options), answer,
        "The mRNA has the sequence of the coding strand with U in place of T. (The polymerase "
        "reads the other strand, the template, to produce it.)",
    ))  # fmt: skip
    protein = "-".join(sequence.three_letter(sequence.translate(EXAMPLE)))
    wrong_frame = "-".join(sequence.three_letter(sequence.translate(EXAMPLE[1:], to_stop=False)))
    options, answer = _place(protein, (wrong_frame, "Met-Lys-Gly-Ser", "Ser-Lys-Gly-Met"), 3)
    q.append(Question(
        "expression",
        f"Translate {rna} with the genetic code, starting at the first letter.",
        options, answer,
        "Read in triplets from the start codon: "
        + ", ".join(f"{sequence.transcribe(c)} = {sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][0]}"
                    for c in sequence.codons(EXAMPLE)) + ".",
    ))  # fmt: skip
    q.append(Question(
        "expression",
        "How many of the 64 codons do not stand for an amino acid?",
        ("1", str(len(sequence.STOP_CODONS)), "4", "20"), 1,
        "The three stop codons "
        + ", ".join(sequence.transcribe(c) for c in sequence.STOP_CODONS)
        + " end translation; the other "
        + str(len(sequence.CODON_TABLE) - len(sequence.STOP_CODONS))
        + " codons name the 20 amino acids.",
    ))  # fmt: skip
    deleted = sequence.apply_mutation(EXAMPLE + "TAA", "deletion", 4)
    effect = sequence.classify_mutation(EXAMPLE + "TAA", deleted).effect
    q.append(Question(
        "mutation",
        "One nucleotide is deleted in the middle of a gene. What is the usual consequence?",
        ("Nothing: the code is degenerate", "One amino acid changes", "A " + effect + ": every "
         "codon after the deletion is read differently", "The gene is no longer replicated"), 2,
        "Codons are read three by three from the start. Removing one letter shifts the reading "
        "frame for the whole rest of the gene.",
    ))  # fmt: skip
    silent = sequence.classify_mutation(
        EXAMPLE + "TAA", sequence.apply_mutation(EXAMPLE + "TAA", "substitution", 5, "C")
    )
    q.append(Question(
        "mutation",
        f"In {EXAMPLE} the codon {silent.codon_before} is changed to {silent.codon_after}. "
        "What happens to the protein?",
        ("It ends early", "Nothing: both codons mean "
         + sequence.AMINO_ACIDS[silent.amino_before][1], "One amino acid changes",
         "The reading frame shifts"), 1,
        "A " + silent.effect + " mutation: the genetic code has several codons for most amino "
        "acids, usually differing in the third letter.",
    ))  # fmt: skip
    return q


def by_topic(topic: str) -> list[Question]:
    return [question for question in questions() if question.topic == topic]
