"""Sequence-level molecular biology: complement, transcription, the genetic code, mutations.

Conventions (the ones of every textbook and of GenBank):

* A DNA sequence is written 5' to 3' and is the **coding strand** of a gene: the mRNA has the same
  letters with U in place of T. The **template strand** is its reverse complement; RNA polymerase
  reads it 3' to 5'.
* Positions are 0-based inside this module; pages and files that show positions to a person use
  1-based positions and say so.

The standard genetic code is built from the compact form published by NCBI (translation table 1,
bases in the order T, C, A, G); a test compares the table with the one shipped in Biopython.
"""

from __future__ import annotations

from dataclasses import dataclass

DNA_LETTERS = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")
HYDROGEN_BONDS = {"A": 2, "T": 2, "G": 3, "C": 3}  # Watson-Crick pair containing this base
BASE_NAMES = {"A": "adenine", "C": "cytosine", "G": "guanine", "T": "thymine", "U": "uracil"}

_NCBI_BASES = "TCAG"
_NCBI_TABLE_1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: dict[str, str] = {
    a + b + c: _NCBI_TABLE_1[16 * i + 4 * j + k]
    for i, a in enumerate(_NCBI_BASES)
    for j, b in enumerate(_NCBI_BASES)
    for k, c in enumerate(_NCBI_BASES)
}
STOP_CODONS = tuple(sorted(c for c, aa in CODON_TABLE.items() if aa == "*"))
START_CODON = "ATG"
# Translation table 11 (bacteria, archaea, phages; the table NCBI gives for phiX174 and lambda) has
# the same codon meanings as table 1 and accepts more start codons.
BACTERIAL_STARTS = ("ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA")

# one-letter code: (three-letter code, name, side-chain class)
AMINO_ACIDS: dict[str, tuple[str, str, str]] = {
    "A": ("Ala", "alanine", "nonpolar"),
    "R": ("Arg", "arginine", "basic"),
    "N": ("Asn", "asparagine", "polar"),
    "D": ("Asp", "aspartate", "acidic"),
    "C": ("Cys", "cysteine", "polar"),
    "Q": ("Gln", "glutamine", "polar"),
    "E": ("Glu", "glutamate", "acidic"),
    "G": ("Gly", "glycine", "nonpolar"),
    "H": ("His", "histidine", "basic"),
    "I": ("Ile", "isoleucine", "nonpolar"),
    "L": ("Leu", "leucine", "nonpolar"),
    "K": ("Lys", "lysine", "basic"),
    "M": ("Met", "methionine", "nonpolar"),
    "F": ("Phe", "phenylalanine", "nonpolar"),
    "P": ("Pro", "proline", "nonpolar"),
    "S": ("Ser", "serine", "polar"),
    "T": ("Thr", "threonine", "polar"),
    "W": ("Trp", "tryptophan", "nonpolar"),
    "Y": ("Tyr", "tyrosine", "polar"),
    "V": ("Val", "valine", "nonpolar"),
    "*": ("Stop", "stop signal", "stop"),
}


def clean_dna(text: str) -> str:
    """Upper-case letters A, C, G, T of a typed sequence (U is read as T, the rest dropped)."""
    return "".join(c for c in text.upper().replace("U", "T") if c in DNA_LETTERS)


def _check(sequence: str) -> str:
    bad = set(sequence) - set(DNA_LETTERS)
    if bad:
        raise ValueError(f"not a DNA sequence: unexpected letters {sorted(bad)}")
    return sequence


def complement(sequence: str) -> str:
    """The partner base of every letter, position by position (so it reads 3' to 5')."""
    return "".join(COMPLEMENT[c] for c in _check(sequence))


def reverse_complement(sequence: str) -> str:
    """The other strand of the double helix, written 5' to 3'."""
    return complement(sequence)[::-1]


def transcribe(coding: str) -> str:
    """The mRNA of a coding strand: the same letters with U in place of T."""
    return _check(coding).replace("T", "U")


def codons(sequence: str) -> list[str]:
    """Consecutive triplets from the first letter; an incomplete last triplet is dropped."""
    return [sequence[i : i + 3] for i in range(0, len(sequence) - 2, 3)]


def translate(coding: str, *, to_stop: bool = True) -> str:
    """Protein (one-letter code) of a coding sequence read from its first letter.

    With ``to_stop`` the protein ends before the first stop codon; else stops are written ``*``.
    """
    protein = []
    for codon in codons(_check(coding)):
        amino = CODON_TABLE[codon]
        if amino == "*" and to_stop:
            break
        protein.append(amino)
    return "".join(protein)


def anticodon(codon: str) -> str:
    """The tRNA anticodon that pairs with an mRNA codon, written 5' to 3' with U."""
    return reverse_complement(codon.replace("U", "T")).replace("T", "U")


def three_letter(protein: str) -> list[str]:
    return [AMINO_ACIDS[a][0] for a in protein]


def base_composition(sequence: str) -> dict[str, int]:
    return {base: sequence.count(base) for base in DNA_LETTERS}


def gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0


def hydrogen_bonds(sequence: str) -> int:
    """Hydrogen bonds that hold this sequence to its complement (2 per A-T, 3 per G-C)."""
    return sum(HYDROGEN_BONDS[c] for c in _check(sequence))


# ---- open reading frames -------------------------------------------------------------------------


@dataclass(frozen=True)
class Orf:
    start: int  # 0-based position of the first letter of the start codon
    end: int  # 0-based position after the last letter of the stop codon
    protein: str

    @property
    def frame(self) -> int:
        return self.start % 3

    @property
    def n_codons(self) -> int:
        return (self.end - self.start) // 3 - 1  # the stop codon is not counted


def find_orfs(
    sequence: str, min_codons: int = 30, starts: tuple[str, ...] = (START_CODON,)
) -> list[Orf]:
    """Open reading frames of the given strand: start codon, in-frame codons, first stop codon.

    Only the longest frame ending at each stop is returned (the first start after the previous
    stop). An open reading frame is a candidate, not a gene: annotation needs other evidence.
    """
    _check(sequence)
    found = []
    for frame in range(3):
        begin: int | None = None
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]
            if begin is None and codon in starts:
                begin = i
            elif begin is not None and CODON_TABLE[codon] == "*":
                if (i - begin) // 3 >= min_codons:
                    protein = "M" + translate(sequence[begin + 3 : i])
                    found.append(Orf(begin, i + 3, protein))
                begin = None
    return sorted(found, key=lambda o: o.start)


# ---- mutations -----------------------------------------------------------------------------------

MUTATION_KINDS = ("substitution", "insertion", "deletion")


def apply_mutation(sequence: str, kind: str, position: int, base: str = "") -> str:
    """The sequence after one point mutation at a 0-based position.

    substitution: the letter at ``position`` becomes ``base``; insertion: ``base`` (one or more
    letters) is inserted before ``position``; deletion: the letter at ``position`` is removed.
    """
    _check(sequence)
    if kind not in MUTATION_KINDS:
        raise ValueError(f"kind must be one of {MUTATION_KINDS}")
    if not 0 <= position < len(sequence):
        raise ValueError("position outside the sequence")
    if kind == "deletion":
        return sequence[:position] + sequence[position + 1 :]
    if not base or set(base) - set(DNA_LETTERS):
        raise ValueError("give the new base(s) as letters A, C, G, T")
    if kind == "insertion":
        return sequence[:position] + base + sequence[position:]
    if len(base) != 1:
        raise ValueError("a substitution replaces one letter by one letter")
    if base == sequence[position]:
        raise ValueError("the new base equals the old one: nothing is mutated")
    return sequence[:position] + base + sequence[position + 1 :]


def substitution_class(old: str, new: str) -> str:
    """transition (purine to purine or pyrimidine to pyrimidine) or transversion."""
    same_family = ({old, new} <= PURINES) or ({old, new} <= PYRIMIDINES)
    return "transition" if same_family else "transversion"


@dataclass(frozen=True)
class MutationEffect:
    effect: str  # silent | missense | nonsense | start lost | stop lost | frameshift | in-frame ...
    codon_index: int | None  # 0-based index of the first changed codon
    codon_before: str
    codon_after: str
    amino_before: str
    amino_after: str
    protein_before: str
    protein_after: str
    explanation: str


def classify_mutation(cds: str, mutated: str) -> MutationEffect:
    """Effect on the protein of changing a coding sequence ``cds`` into ``mutated``.

    ``cds`` must start with a start codon and end with a stop codon. The classification follows
    the usual textbook definitions: *silent* (same amino acid), *missense* (another amino acid),
    *nonsense* (an amino-acid codon becomes a stop), *frameshift* (length changes by a number of
    letters not divisible by 3, so every later codon is read differently), *in-frame* insertion
    or deletion (length changes by a multiple of 3).
    """
    before = translate(cds)
    after = translate(mutated)
    change = len(mutated) - len(cds)
    old_codons, new_codons = codons(cds), codons(mutated)
    index = next(
        (k for k, (a, b) in enumerate(zip(old_codons, new_codons, strict=False)) if a != b), None
    )
    c_old = old_codons[index] if index is not None else ""
    c_new = new_codons[index] if index is not None else ""
    a_old = CODON_TABLE.get(c_old, "")
    a_new = CODON_TABLE.get(c_new, "")

    def result(effect: str, explanation: str) -> MutationEffect:
        return MutationEffect(effect, index, c_old, c_new, a_old, a_new, before, after, explanation)

    if change % 3 != 0:
        return result(
            "frameshift",
            f"The length changes by {change:+d} letter(s), not a multiple of 3: from codon "
            f"{(index or 0) + 1} on, every triplet is read in a new frame. The protein has "
            f"{len(after)} amino acids instead of {len(before)}.",
        )
    if change != 0:
        kind = "insertion" if change > 0 else "deletion"
        return result(
            f"in-frame {kind}",
            f"The length changes by {change:+d} letters, a multiple of 3: {abs(change) // 3} amino "
            f"acid(s) are {'added' if change > 0 else 'lost'} and the reading frame is kept.",
        )
    if index is None:
        return result("no change", "The two sequences are identical.")
    name_old, name_new = AMINO_ACIDS[a_old][0], AMINO_ACIDS[a_new][0]
    where = f"Codon {index + 1}: {c_old} ({name_old}) becomes {c_new} ({name_new})."
    if index == 0 and c_new not in BACTERIAL_STARTS:
        return result("start lost", where + " The start codon is destroyed: no protein is begun.")
    if a_old == a_new:
        return result(
            "silent",
            where + " Both codons mean the same amino acid (the code is degenerate), so the "
            "protein is unchanged.",
        )
    if a_new == "*":
        return result(
            "nonsense",
            where + f" A stop codon appears early: the protein is cut to {len(after)} of "
            f"{len(before)} amino acids.",
        )
    if a_old == "*":
        return result(
            "stop lost", where + " The stop codon is destroyed: the ribosome reads past the end."
        )
    class_old, class_new = AMINO_ACIDS[a_old][2], AMINO_ACIDS[a_new][2]
    nature = (
        f"both are {class_old}, a conservative change"
        if class_old == class_new
        else f"{class_old} becomes {class_new}, a non-conservative change"
    )
    return result("missense", where + f" One amino acid of the protein changes ({nature}).")
