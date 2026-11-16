"""Scenarios of the three animated processes: replication, transcription, translation.

A scenario is the data of an animation: the real sequence it uses, everything derived from it
(complement, RNA, codons, amino acids) and the ordered steps with their explanations. The player
in ``genomeasm/appkit/web/process.js`` only moves these pieces; it decides no biology.

Each step has a ``to`` dictionary: the state of the scene at the end of the step. The player
interpolates between the states of consecutive steps, so any moment of the animation can be
reached by scrubbing and every frame is a function of time only.

The drawings are schematic and say so on screen: the helix is drawn flat, proteins are plain
shapes, and lengths are shortened (an Okazaki fragment here has 8 nucleotides).
"""

from __future__ import annotations

from typing import Any

from genomeasm.molbio import library, sequence

REPLICATION_LENGTH = 24
OKAZAKI_LENGTH = 8
PRIMER_LENGTH = 3
LAMBDA_WINDOW_START = 227  # 1-based; a stretch of lambda with all four letters well mixed
TRANSCRIPT_LENGTH = 24
TRANSCRIPTION_FLANK = 3


Window = dict[str, tuple[float, float]]


def _step(
    title: str, text: str, seconds: float, win: Window | None = None, **to: float
) -> dict[str, Any]:
    """One step. ``win`` gives, for some keys, the part of the step (0 to 1) in which that
    value changes; the other values change over the whole step."""
    return {"title": title, "text": text, "dur": seconds, "to": to, "win": win or {}}


def replication_scenario(genome: str, start: int = LAMBDA_WINDOW_START,
                         source: str = "lambda, NC_001416.1") -> dict[str, Any]:  # fmt: skip
    """A replication fork passing over ``REPLICATION_LENGTH`` base pairs of a real sequence.

    The top strand is written 5' to 3' from left to right and the fork moves to the right. The
    new strand on the bottom template therefore grows towards the fork (leading strand) and the
    new strand on the top template away from it, in fragments (lagging strand).
    """
    top = genome[start - 1 : start - 1 + REPLICATION_LENGTH]
    if len(top) != REPLICATION_LENGTH:
        raise ValueError("the window does not fit in the genome")
    bottom = sequence.complement(top)
    n = REPLICATION_LENGTH
    fragments = [[a, a + OKAZAKI_LENGTH - 1] for a in range(0, n, OKAZAKI_LENGTH)]
    o = float(OKAZAKI_LENGTH)
    steps = [
        _step(
            "The parental double helix",
            "Two <b>antiparallel</b> strands: the upper one runs 5′→3′ from left to right, the "
            "lower one 3′→5′. Because A pairs only with T and G only with C, <b>each strand "
            "determines the other</b>, so each can serve as the template for a new partner.",
            5, fork=0, lead=0, f0=0, f1=0, f2=0, swap=0, ligate=0, apart=0,
        ),
        _step(
            "Helicase opens the helix",
            "The enzyme <b>helicase</b> breaks the hydrogen bonds between the base pairs and "
            "unwinds the helix. The opening point is the <b>replication fork</b>. Single-strand "
            "binding proteins keep the separated strands from pairing again.",
            5, fork=10, lead=0, f0=0, f1=0, f2=0, swap=0, ligate=0, apart=0,
        ),
        _step(
            "Primase makes RNA primers",
            "DNA polymerase cannot begin a strand, it can only lengthen one. <b>Primase</b> "
            "therefore builds a short <b>RNA primer</b> on each template (RNA uses U where DNA "
            "has T). The primer offers the free 3′ end the polymerase needs.",
            5, fork=10, lead=PRIMER_LENGTH, f0=PRIMER_LENGTH, f1=0, f2=0, swap=0, ligate=0, apart=0,
        ),
        _step(
            "DNA polymerase: always 5′ → 3′",
            "<b>DNA polymerase</b> adds each new nucleotide to the 3′ end, so a new strand always "
            "grows <b>5′→3′</b>. On the lower template that direction follows the fork: this "
            "<b>leading strand</b> is made in one piece.",
            6, fork=10, lead=8, f0=o, f1=0, f2=0, swap=0, ligate=0, apart=0,
        ),
        _step(
            "The lagging strand is made in pieces",
            "On the upper template 5′→3′ points <b>away from the fork</b>. This <b>lagging "
            "strand</b> is made backwards in short pieces, the <b>Okazaki fragments</b>; each "
            "new piece starts with its own primer as soon as the fork has opened enough template.",
            8, {"fork": (0.0, 0.5), "lead": (0.0, 0.5), "f1": (0.5, 1.0)},
            fork=18, lead=16, f0=o, f1=o, f2=0, swap=0, ligate=0, apart=0,
        ),
        _step(
            "The fork moves on",
            "The leading strand keeps growing without a break while a third Okazaki fragment is "
            "made. One fork thus works with two polymerases, one per template.",
            8, {"fork": (0.0, 0.5), "lead": (0.0, 0.5), "f2": (0.5, 1.0)},
            fork=float(n + 4), lead=float(n), f0=o, f1=o, f2=o, swap=0, ligate=0, apart=0,
        ),
        _step(
            "Primers are replaced, ligase joins the pieces",
            "The RNA primers are removed and the gaps filled with DNA (in bacteria by DNA "
            "polymerase I). <b>DNA ligase</b> then closes the remaining nicks: the lagging strand "
            "becomes one continuous strand.",
            6, {"swap": (0.0, 0.5), "ligate": (0.5, 1.0)},
            fork=float(n + 4), lead=float(n), f0=o, f1=o, f2=o, swap=1, ligate=1, apart=0,
        ),
        _step(
            "Two identical helices: semiconservative replication",
            "Each daughter helix has <b>one old strand and one new strand</b>. This "
            "<b>semiconservative</b> mechanism was demonstrated by Meselson and Stahl in 1958. "
            "Both copies carry the sequence of the original.",
            6, fork=float(n + 4), lead=float(n), f0=o, f1=o, f2=o, swap=1, ligate=1, apart=1,
        ),
    ]  # fmt: skip
    return {
        "kind": "replication",
        "title": "DNA replication",
        "note": f"Schematic. Sequence: positions {start:,} to {start + n - 1:,} of {source}. The "
        "helix is drawn flat; real Okazaki fragments are far longer than the 8 nucleotides shown.",
        "data": {
            "top": top,
            "bottom": bottom,
            "fragments": fragments,
            "primer": PRIMER_LENGTH,
            "start": start,
        },
        "steps": steps,
    }


def transcription_scenario(genome: str, gene_name: str = "J") -> dict[str, Any]:
    """RNA polymerase transcribing the first ``TRANSCRIPT_LENGTH`` letters of a phiX174 gene."""
    gene = library.gene(gene_name)
    first = gene.parts[0][0]  # 1-based position of the A of ATG
    flank, n = TRANSCRIPTION_FLANK, TRANSCRIPT_LENGTH
    coding = genome[first - 1 - flank : first - 1 + n + flank]
    template = sequence.complement(coding)
    rna = sequence.transcribe(coding[flank : flank + n])
    length = float(n)
    steps = [
        _step(
            "A gene on the double helix",
            "A gene is a stretch of DNA that is copied into RNA. The <b>coding strand</b> (top) "
            "has the sequence the RNA will have; the <b>template strand</b> (bottom) is the one "
            "that is actually read.",
            5, enzyme=0, bubble=0, pos=0, release=0, compare=0,
        ),
        _step(
            "Initiation: RNA polymerase opens the helix",
            "<b>RNA polymerase</b> binds at a <b>promoter</b> and separates the two strands over "
            "a short stretch, the <b>transcription bubble</b>. Unlike DNA polymerase it needs "
            "no primer.",
            5, enzyme=1, bubble=1, pos=0, release=0, compare=0,
        ),
        _step(
            "Elongation: reading 3′ → 5′, writing 5′ → 3′",
            "The polymerase reads the template strand 3′→5′ and joins ribonucleotides 5′→3′. "
            "The pairing rule is the same as in DNA, except that RNA carries <b>uracil (U)</b> "
            "instead of thymine: an A of the template pairs with U.",
            7, enzyme=1, bubble=1, pos=8, release=0, compare=0,
        ),
        _step(
            "The bubble travels with the enzyme",
            "Only a few nucleotides of the new RNA stay paired with the template. Behind the "
            "polymerase the DNA closes again and the RNA leaves as a <b>single strand</b>.",
            7, enzyme=1, bubble=1, pos=18, release=0, compare=0,
        ),
        _step(
            "Termination: the RNA is released",
            "At a <b>terminator</b> the polymerase lets go of the DNA and of the finished RNA; "
            "the double helix is intact again and can be transcribed many times.",
            6, enzyme=0, bubble=0, pos=length, release=1, compare=0,
        ),
        _step(
            "The transcript equals the coding strand, with U for T",
            "Compare the two rows: the RNA repeats the coding strand letter by letter, <b>U in "
            "place of T</b>. In bacteria this messenger RNA is translated at once; in eukaryotes "
            "it is first capped, spliced and given a poly-A tail.",
            6, enzyme=0, bubble=0, pos=length, release=1, compare=1,
        ),
    ]  # fmt: skip
    return {
        "kind": "transcription",
        "title": "Transcription: DNA to RNA",
        "note": f"Schematic. Sequence: the first {n} letters of gene {gene_name} of phiX174 "
        f"(positions {first:,} to {first + n - 1:,} of NC_001422.1). In the phage this gene lies "
        "inside a longer transcript; the promoter is drawn next to it for teaching.",
        "data": {
            "coding": coding,
            "template": template,
            "rna": rna,
            "flank": flank,
            "first": first,
            "gene": gene_name,
        },  # fmt: skip
        "steps": steps,
    }


def translation_scenario(genome: str, gene_name: str = "J") -> dict[str, Any]:
    """A ribosome translating a whole (short) gene of phiX174 from its real sequence."""
    gene = library.gene(gene_name)
    cds = gene.sequence(genome)
    triplets = sequence.codons(cds)
    protein = sequence.translate(cds)
    if len(triplets) != len(protein) + 1 or sequence.CODON_TABLE[triplets[-1]] != "*":
        raise ValueError(f"gene {gene_name} is not a clean coding sequence")
    codon_rows = [
        {
            "codon": sequence.transcribe(c),
            "anticodon": sequence.complement(c).replace("T", "U"),  # 3' to 5', under the codon
            "aa": sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][0],
            "one": sequence.CODON_TABLE[c],
            "cls": sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][2],
        }
        for c in triplets
    ]
    last = float(len(protein) - 1)  # index of the last amino-acid codon
    steps = [
        _step(
            "The messenger RNA is read in codons",
            "Three nucleotides, one <b>codon</b>, stand for one amino acid. The 64 codons are "
            "the <b>genetic code</b>: 61 mean an amino acid, 3 mean stop. The <b>start codon "
            "AUG</b> fixes the reading frame.",
            5, ribo=0, cycle=0, rf=0, free=0,
        ),
        _step(
            "Initiation at the start codon",
            "The small ribosomal subunit binds the mRNA, the initiator tRNA carrying "
            "<b>methionine</b> pairs with AUG, then the large subunit joins. The ribosome has "
            "three tRNA sites: <b>A</b> (arrival), <b>P</b> (peptide) and <b>E</b> (exit).",
            5, ribo=1, cycle=0, rf=0, free=0,
        ),
        _step(
            "A tRNA brings the next amino acid",
            "A <b>transfer RNA</b> carries one amino acid and has an <b>anticodon</b> of three "
            "nucleotides. Only the tRNA whose anticodon pairs with the codon in the A site is "
            "accepted: this is how a nucleotide sequence is translated into an amino-acid sequence.",
            6, ribo=1, cycle=0.4, rf=0, free=0,
        ),
        _step(
            "A peptide bond forms",
            "The growing chain is passed from the tRNA in the P site onto the new amino acid in "
            "the A site. The catalyst is ribosomal RNA, not protein: <b>the ribosome is a "
            "ribozyme</b> (Nissen et al. 2000).",
            5, ribo=1, cycle=0.6, rf=0, free=0,
        ),
        _step(
            "Translocation: three nucleotides forward",
            "The ribosome moves one codon along the mRNA in the 5′→3′ direction. The emptied "
            "tRNA leaves through the E site and the A site is free for the next tRNA.",
            5, ribo=1, cycle=1, rf=0, free=0,
        ),
        _step(
            "Elongation, codon after codon",
            "The three moves repeat for every codon. The protein grows from its amino end "
            "(N-terminus) to its carboxyl end (C-terminus), one amino acid per three nucleotides.",
            14, ribo=1, cycle=last, rf=0, free=0,
        ),
        _step(
            "Termination at a stop codon",
            "No tRNA matches the stop codons UAA, UAG and UGA. A protein <b>release factor</b> "
            "enters the A site instead, and the finished chain is cut loose from the last tRNA.",
            6, ribo=1, cycle=last, rf=1, free=0,
        ),
        _step(
            f"The finished protein: {len(protein)} amino acids",
            f"The ribosome separates into its subunits. The product is protein {gene_name} of "
            f"phiX174 ({gene.product}), translated here codon by codon from the real gene.",
            6, ribo=0, cycle=last, rf=1, free=1,
        ),
    ]  # fmt: skip
    return {
        "kind": "translation",
        "title": "Translation: RNA to protein",
        "note": f"Schematic. Sequence: gene {gene_name} of phiX174, positions "
        f"{gene.parts[0][0]:,} to {gene.parts[-1][1]:,} of NC_001422.1 ({len(cds)} nucleotides, "
        f"{len(protein)} amino acids, NCBI protein {gene.protein_id}).",
        "data": {"codons": codon_rows, "protein": protein, "gene": gene_name},
        "steps": steps,
    }


def scenarios(lambda_genome: str, phix_genome: str) -> dict[str, dict[str, Any]]:
    return {
        "replication": replication_scenario(lambda_genome),
        "transcription": transcription_scenario(phix_genome),
        "translation": translation_scenario(phix_genome),
    }
