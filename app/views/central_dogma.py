"""From DNA to protein: replication, transcription and translation as step-by-step animations."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import biology, studio  # noqa: E402
from genomeasm.molbio import diagrams, library, sequence  # noqa: E402

st.title("From DNA to protein", icon=":material/conversion_path:")
st.markdown(
    "DNA is **replicated** so that every cell inherits a copy, **transcribed** into RNA, and the RNA is "
    "**translated** into protein. The three animations below use real sequences of the two genomes of "
    "this project; play them whole or step by step."
)
section = ui.section_picker(
    "Section",
    [
        "Overview",
        "Replication",
        "Transcription",
        "Genetic code",
        "Translation",
        "Try it",
        "Check yourself",
    ],
    "dogma-section",
)
scenarios = biology.scenarios()


def animation(kind: str) -> None:
    ui.studio(
        studio.process_html(scenarios[kind], height=700),
        720,
        f"Animation: {scenarios[kind]['title']}",
    )
    st.caption("Keys while the animation has the focus: space = play or pause, arrows = one step.")


if section == "Overview":
    st.subheader("The central dogma of molecular biology")
    ui.diagram(diagrams.central_dogma_svg(), "The central dogma: DNA to RNA to protein")
    st.markdown(
        "- **Replication**: DNA → DNA. Before a cell divides, the whole genome is copied.\n"
        "- **Transcription**: DNA → RNA. A gene is copied into a messenger RNA.\n"
        "- **Translation**: RNA → protein. The ribosome reads the messenger RNA in codons of three "
        "nucleotides and joins the matching amino acids.\n"
        "- **RNA** differs from DNA in three ways: its sugar is ribose, it uses **uracil (U)** instead of "
        "thymine, and it is usually single-stranded."
    )
    ui.takeaway(library.fact_text("central_dogma"))
    st.caption(
        "The dashed arrows are the special transfers named by Crick: retroviruses copy RNA into DNA, and "
        "some RNA viruses copy RNA directly."
    )
    ui.sources("crick1970", "alberts2022")

elif section == "Replication":
    st.subheader("Replication: each strand is the template for a new one")
    ui.key_facts(
        [
            ui.fact_item("semiconservative"),
            ui.fact_item("fork_speed"),
            ui.fact_item("fidelity"),
            ui.fact_item("okazaki"),
        ]
    )
    animation("replication")
    with st.expander("The enzymes at the replication fork", icon=":material/science:"):
        ui.table(
            [
                {
                    "Protein": "Helicase",
                    "What it does": "Breaks the hydrogen bonds and unwinds the helix at the fork",
                },
                {
                    "Protein": "Single-strand binding protein",
                    "What it does": "Keeps the separated strands apart",
                },
                {
                    "Protein": "Primase",
                    "What it does": "Makes a short RNA primer, the starting point for DNA polymerase",
                },
                {
                    "Protein": "DNA polymerase",
                    "What it does": "Adds nucleotides to the 3′ end of the new strand (5′ → 3′) and proofreads",
                },
                {
                    "Protein": "DNA ligase",
                    "What it does": "Joins the Okazaki fragments into one strand",
                },
            ]
        )
        st.caption("Names as used for bacteria; eukaryotes use more proteins with the same roles.")
    ui.takeaway(
        "DNA polymerase only extends 3′ ends, and the two templates are antiparallel. So one new strand is "
        "made continuously (leading strand) and the other in fragments (lagging strand)."
    )
    ui.sources("meselson1958", "okazaki1968", "alberts2022", "refseq_lambda")

elif section == "Transcription":
    st.subheader("Transcription: a gene is copied into RNA")
    data = scenarios["transcription"]["data"]
    st.markdown(
        "RNA polymerase opens the helix, reads the **template strand** 3′ → 5′ and builds an RNA 5′ → 3′. "
        "The RNA therefore has the sequence of the other strand, the **coding strand**, with U for T."
    )
    animation("transcription")
    flank = data["flank"]
    coding = data["coding"][flank : flank + len(data["rna"])]
    st.code(
        f"coding strand    5′ {coding} 3′\n"
        f"template strand  3′ {sequence.complement(coding)} 5′\n"
        f"mRNA             5′ {data['rna']} 3′",
        language=None,
    )
    ui.takeaway(
        "Only one strand of a gene is read. The messenger RNA is a mobile, single-stranded copy of the "
        "coding strand."
    )
    ui.sources("refseq_phix174", "alberts2022")

elif section == "Genetic code":
    st.subheader("The genetic code: three nucleotides name one amino acid")
    n_stop = len(sequence.STOP_CODONS)
    ui.key_facts(
        [
            ("Codons", str(len(sequence.CODON_TABLE)), "4 × 4 × 4 triplets of four letters"),
            (
                "Codons for amino acids",
                str(len(sequence.CODON_TABLE) - n_stop),
                "Counted from the code table of this project",
            ),
            (
                "Stop codons",
                ", ".join(sequence.transcribe(c) for c in sequence.STOP_CODONS),
                "They end translation",
            ),
            ("Start codon", "AUG", "Methionine; also sets the reading frame"),
        ]
    )
    ui.studio(
        studio.fragment_html(diagrams.codon_table_html(), diagrams.CODON_TABLE_CSS),
        560,
        "The standard genetic code as a table of 64 codons",
    )
    st.caption(
        "First letter: row group; second letter: column; third letter: right-hand column. Colours: class of "
        "the amino-acid side chain (grey nonpolar, blue polar, orange acidic, green basic, red stop)."
    )
    st.markdown(
        "- The code is **degenerate**: most amino acids have several codons, usually differing in the "
        "**third** letter. This is why many single-letter mutations change nothing.\n"
        f"- {library.fact_text('triplet')}\n"
        "- The code is nearly the same in all organisms. The two viruses of this project use NCBI "
        "translation table 11 (bacteria), which has the same codon meanings and allows more start codons."
    )
    ui.sources("crick1961", "nirenberg1961", "ncbi_gencode")

elif section == "Translation":
    st.subheader("Translation: the ribosome builds a protein")
    gene = library.gene("J")
    protein = sequence.translate(biology.gene_sequence("J"))
    ui.key_facts(
        [
            (
                "Gene shown",
                f"J of phiX174, {gene.length} nt",
                f"NCBI protein {gene.protein_id}, positions {gene.parts[0][0]} to {gene.parts[-1][1]}",
            ),
            (
                "Protein",
                f"{len(protein)} amino acids",
                "Translated here from the sequence in data/raw",
            ),
            ui.fact_item("ribozyme"),
        ]
    )
    animation("translation")
    st.code(f"protein J: {protein}", language=None)
    ui.takeaway(
        "A transfer RNA couples one amino acid to one anticodon. Codon by codon, the ribosome turns a "
        "sequence of nucleotides into a sequence of amino acids."
    )
    ui.sources("nissen2000", "refseq_phix174", "alberts2022")

elif section == "Try it":
    st.subheader("Translate a sequence yourself")
    genes = library.phix_genes()
    names = [g.name for g in genes]
    c1, c2 = st.columns([1, 3])
    picked = c1.selectbox(
        "Gene of phiX174",
        names,
        index=names.index("J"),
        format_func=lambda n: f"{n}: {library.gene(n).product}",
        key="try-gene",
    )
    typed = c2.text_input(
        "or type your own DNA (coding strand, 5′ → 3′)",
        key="try-dna",
        placeholder="for example ATGGCCAAGTAA",
    )
    dna = sequence.clean_dna(typed) if typed.strip() else biology.gene_sequence(picked)
    if len(dna) < 3:
        st.warning("Type at least three letters A, C, G or T.")
        st.stop()
    protein = sequence.translate(dna)
    stops = sequence.translate(dna, to_stop=False).count("*")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Nucleotides", f"{len(dna):,}", border=True)
    m2.metric(
        "Codons read",
        f"{len(protein) + (1 if stops else 0):,}",
        border=True,
        help="Up to and including the first stop codon",
    )
    m3.metric("Amino acids", f"{len(protein):,}", border=True)
    m4.metric("GC content", f"{100 * sequence.gc_fraction(dna):.1f}%", border=True)
    if dna[:3] not in sequence.BACTERIAL_STARTS:
        st.warning(
            f"The sequence does not begin with a start codon (it begins with {dna[:3]}); it is translated from its first letter anyway."
        )
    if not stops:
        st.info(
            "No stop codon in this reading frame: the ribosome would run to the end of the sequence."
        )
    shown = dna[:90]
    st.code(
        f"DNA (coding)  5′ {shown}{'…' if len(dna) > 90 else ''}\n"
        f"template      3′ {sequence.complement(shown)}{'…' if len(dna) > 90 else ''}\n"
        f"mRNA          5′ {sequence.transcribe(shown)}{'…' if len(dna) > 90 else ''}\n"
        f"protein          {protein}",
        language=None,
    )
    rows = [
        {
            "codon": k + 1,
            "DNA": c,
            "mRNA": sequence.transcribe(c),
            "anticodon (5′→3′)": sequence.anticodon(c),
            "amino acid": sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][1],
            "code": sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][0],
            "class": sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][2],
        }
        for k, c in enumerate(sequence.codons(dna)[: len(protein) + 1])
    ]
    ui.table(rows, height=320)
    st.caption(
        "Computed live with the code table of this project (compared with Biopython by a test)."
    )

else:
    st.subheader("Check yourself: replication and expression")
    ui.quiz("replication")
    ui.quiz("expression")
