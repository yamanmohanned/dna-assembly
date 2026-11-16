"""Genes, genomes and mutations: genome sizes, packaging, the gene map of phiX174, a mutation lab."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import biology  # noqa: E402
from genomeasm.molbio import diagrams, library, sequence  # noqa: E402

st.title("Genes, genomes and mutations", icon=":material/account_tree:")
st.markdown(
    "A **gene** is a stretch of DNA that is copied into RNA. The **genome** is the complete DNA of an "
    "organism, carried on one or more **chromosomes**. A **mutation** is a heritable change of its sequence."
)
section = ui.section_picker(
    "Section",
    [
        "Genomes",
        "Packaging",
        "Gene map of phiX174",
        "Mutation lab",
        "Biotechnology",
        "Sequencing",
        "Check yourself",
    ],
    "genome-section",
)
facts = library.facts()
phix = biology.genome(biology.PHIX)

if section == "Genomes":
    st.subheader("Genomes differ in size by a factor of a million")
    ui.key_facts(
        [
            ui.fact_item("human"),
            ui.fact_item("human_genes"),
            ui.fact_item("chromosomes_human"),
            ui.fact_item("ecoli"),
        ]
    )
    rows = biology.genome_rows()
    ui.diagram(
        diagrams.genome_sizes_svg(rows),
        "Genome sizes on a logarithmic scale",
        "Dark bars: the two genomes of this project (counted from data/raw). The other values are quoted "
        "from the papers that reported the sequences.",
    )
    ui.table(
        [{"Genome": name, "Base pairs": f"{size:,}", "Source": note} for name, size, note in rows]
    )
    st.markdown(
        f"A human genome is about **{rows[-1][1] // rows[0][1]:,} times** longer than that of phiX174. Size is "
        "not complexity in a simple sense: most of a human genome does not code for protein, while the "
        "small viral genomes are packed with genes (next sections)."
    )
    ui.takeaway(
        "Small viral genomes are ideal for studying assembly: they are complete, public and small enough "
        "to repeat every experiment on a laptop."
    )
    ui.sources("fleischmann1995", "blattner1997", "goffeau1996", "nurk2022", "ihgsc2004")

elif section == "Packaging":
    st.subheader("How a long DNA molecule fits into a cell")
    ui.diagram(
        diagrams.packaging_svg(),
        "Levels of DNA packaging in a eukaryotic cell",
        "Schematic, not to scale.",
    )
    ui.key_facts(
        [
            ui.fact_item("nucleosome"),
            ui.fact_item("chromatin_fibre"),
            ui.fact_item("chromosomes_human"),
        ]
    )
    st.markdown(
        f"- **Nucleosomes.** {facts['nucleosome']['statement']} Seen in the electron microscope they look "
        "like beads on a string.\n"
        f"- **Higher levels.** {facts['chromatin_fibre']['statement']} This is a point where textbooks "
        "changed: the neat 30-nm fibre is not what is found inside cells.\n"
        "- **Chromosomes.** The most compact state is reached during cell division, when each replicated "
        "chromosome consists of two sister chromatids joined at the centromere.\n"
        "- **Bacteria and viruses** have no histones. A bacterial chromosome is compacted by other proteins "
        "and supercoiling; a phage genome is packed tightly inside its protein shell."
    )
    ui.takeaway(
        "Packaging is not only storage: how tightly a region is packed helps decide whether its genes are read."
    )
    ui.sources("kornberg1974", "luger1997", "ou2017", "maeshima2010")

elif section == "Gene map of phiX174":
    st.subheader("A complete genome: eleven genes in 5,386 nucleotides")
    genes = library.phix_genes()
    summary = biology.overlap_summary()
    ui.key_facts(
        [
            ("Protein-coding genes", str(summary["n_genes"]), "NCBI feature table of NC_001422.1"),
            ("Genome length", f"{summary['genome_length']:,} nt", "Counted from data/raw"),
            (
                "Sum of the gene lengths",
                f"{summary['sum_of_gene_lengths']:,} nt",
                "More than the genome: genes overlap",
            ),
            (
                "Letters inside a gene",
                f"{100 * summary['coding_letters'] / summary['genome_length']:.0f}%",
                "Share of the genome covered by at least one gene",
            ),
        ]
    )
    ui.diagram(
        diagrams.gene_map_svg(genes, len(phix), library.phix_annotation()["transcripts"]),
        "Gene map of bacteriophage phiX174",
        "One row per reading frame. Genes that share DNA lie on different rows or sub-rows. ↩ marks genes "
        "that run over the end of the record: the genome is a circle.",
    )
    ui.table(
        [
            {
                "Gene": g.name,
                "Position": " + ".join(f"{a:,}–{b:,}" for a, b in g.parts),
                "Length (nt)": g.length,
                "Protein (amino acids)": len(sequence.translate(g.sequence(phix))),
                "Frame": g.frame() + 1,
                "Product (NCBI)": g.product,
                "Protein ID": g.protein_id,
            }
            for g in genes
        ]
    )
    pairs = ", ".join(f"{a} and {b} ({n:,} nt)" for a, b, n in summary["pairs"] if n > 4)
    st.markdown(
        f"**Overlapping genes.** {facts['phix_overlap']['statement']} Counted from the annotation, the "
        f"larger overlaps are: {pairs}. Gene A* is a shorter protein made from an internal start inside gene A, "
        "in the same frame; genes B, K and E are read in a different frame from the gene they lie in."
    )
    ui.takeaway(library.fact_text("phix_first"))
    st.caption(
        "Every gene in this table is checked by a test against the stored sequence: start codon, whole "
        "codons, exactly one stop codon at the end."
    )
    ui.sources("sanger1977phix", "barrell1976", "refseq_phix174")

elif section == "Mutation lab":
    st.subheader("Change one letter of a real gene and see what happens to the protein")
    genes = library.phix_genes()
    names = [g.name for g in genes]
    c1, c2, c3, c4 = st.columns([2, 2, 2, 2], vertical_alignment="bottom")
    name = c1.selectbox(
        "Gene of phiX174", names, index=names.index("J"),
        format_func=lambda n: f"{n}: {library.gene(n).product}", key="mut-gene",
    )  # fmt: skip
    cds = biology.gene_sequence(name)
    kind = c2.segmented_control(
        "Kind of mutation",
        list(sequence.MUTATION_KINDS),
        default="substitution",
        key="mut-kind",
        required=True,
    )
    position = int(c3.number_input("Letter of the gene", 1, len(cds), 7, 1, key="mut-position"))
    position = min(position, len(cds))
    old = cds[position - 1]
    base = ""
    if kind != "deletion":
        choices = [b for b in "ACGT" if b != old] if kind == "substitution" else list("ACGT")
        base = str(
            c4.segmented_control(
                "New base" if kind == "substitution" else "Inserted base", choices,
                default=choices[0], key=f"mut-base-{kind}-{old}", required=True,
            )
        )  # fmt: skip
    mutated = sequence.apply_mutation(cds, str(kind), position - 1, base)
    effect = sequence.classify_mutation(cds, mutated)
    tone = {"silent": "green", "no change": "green", "missense": "orange"}.get(effect.effect, "red")
    st.badge(effect.effect.upper(), color=tone, icon=":material/genetics:")
    st.markdown(f"**{effect.explanation}**")
    if kind == "substitution":
        st.caption(
            f"{old} → {base} is a {sequence.substitution_class(old, base)} (purine ↔ purine or pyrimidine ↔ pyrimidine is a transition)."
        )
    codon_index = (position - 1) // 3
    start = max(0, codon_index - 4) * 3
    window = 36

    def spaced(text: str) -> str:
        return " ".join(text[i : i + 3] for i in range(0, len(text), 3))

    def amino_line(text: str) -> str:
        return " ".join(
            f"{sequence.AMINO_ACIDS[sequence.CODON_TABLE[c]][0]}" for c in sequence.codons(text)
        )

    before, after = cds[start : start + window], mutated[start : start + window]
    offset = position - 1 - start  # the changed letter inside the window; codons are spaced
    marker = " " * (offset + offset // 3) + "^"
    st.code(
        f"codon number    {start // 3 + 1} …\n"
        f"before  DNA     {spaced(before)}\n"
        f"                {marker}\n"
        f"        protein {amino_line(before)}\n"
        f"after   DNA     {spaced(after)}\n"
        f"        protein {amino_line(after)}",
        language=None,
    )
    m1, m2, m3 = st.columns(3)
    m1.metric("Protein before", f"{len(effect.protein_before)} amino acids", border=True)
    m2.metric(
        "Protein after", f"{len(effect.protein_after)} amino acids", border=True,
        delta=len(effect.protein_after) - len(effect.protein_before), delta_color="off",
    )  # fmt: skip
    same = sum(a == b for a, b in zip(effect.protein_before, effect.protein_after, strict=False))
    m3.metric(
        "Amino acids unchanged",
        f"{same} of {len(effect.protein_before)}",
        border=True,
        help="Compared position by position",
    )
    with st.expander("The whole protein, before and after"):
        st.code(f"before: {effect.protein_before}\nafter:  {effect.protein_after}", language=None)
    with st.expander("Worked examples on gene J", icon=":material/list:"):
        ui.table(
            [
                {
                    "Change": r["change"],
                    "Codon": r["codon"],
                    "Effect": r["effect"],
                    "What happens": r["explanation"],
                }
                for r in biology.mutation_rows()
            ]
        )
    st.markdown(
        f"**In medicine.** {library.fact_text('sickle')} "
        "**Do not confuse** mutations with the *read errors* of the Simulator page: a read error is a mistake "
        "of the sequencing machine in one read, not a change in the organism's DNA."
    )
    ui.takeaway(
        "The effect of a mutation is decided by the genetic code: a substitution changes at most one amino "
        "acid, an insertion or deletion that is not a multiple of three changes everything after it."
    )
    ui.sources("ingram1957", "refseq_phix174", "ncbi_gencode")

elif section == "Biotechnology":
    st.subheader("Reading, copying and editing DNA")
    st.markdown("Modern biotechnology uses the enzymes that cells themselves use on DNA.")
    for key, title, icon in (
        ("sanger_method", "Sequencing: reading DNA", ":material/menu_book:"),
        ("pcr", "PCR: copying a chosen stretch", ":material/content_copy:"),
        ("crispr", "CRISPR-Cas9: cutting at a chosen site", ":material/content_cut:"),
    ):
        with st.container(border=True):
            st.markdown(f"#### {icon} {title}")
            st.markdown(f"**{facts[key]['value']}.** {facts[key]['statement']}")
            st.caption(
                f"{library.cite(str(facts[key]['source']))} · {library.CHECK_LABELS[str(facts[key]['check'])]}"
            )
    st.markdown(
        "All three rest on **base pairing**: a primer finds its target in PCR, a guide RNA finds its target "
        "for Cas9, and in sequencing a polymerase copies a template letter by letter."
    )
    ui.sources("sanger1977seq", "saiki1988", "jinek2012")

elif section == "Sequencing":
    st.subheader("From sequencing to the assembly problem")
    left, right = st.columns([2, 3])
    letters, first = biology.sanger_window()
    with left:
        ui.diagram(
            diagrams.sanger_gel_svg(letters, first),
            "Idealised chain-termination sequencing gel",
            "Idealised gel computed from the real sequence at the start of gene J.",
        )
    with right:
        st.markdown(
            f"**Chain-termination sequencing.** {facts['sanger_method']['statement']} Each of the four "
            "lanes contains copies that end with one known letter; the shortest copies run furthest. Reading "
            "the bands from the bottom gives the sequence.\n\n"
            "**The limit.** One reaction reads only a short stretch. A genome is therefore broken into many "
            "random pieces that are sequenced separately: **shotgun sequencing**. "
            f"{facts['hinfluenzae']['statement']}"
        )
    shot = biology.shotgun_example()
    ui.diagram(
        diagrams.shotgun_svg(shot["starts"], shot["lengths"], shot["length"]),
        "Shotgun reads at their true positions and their coverage",
        f"A real simulation: {shot['n']} reads of {shot['read_length']} letters from the first {shot['length']} "
        f"letters of phiX174 at coverage {shot['coverage']:.0f}×. The sequencer does not report the positions: "
        "finding them again is the assembly problem.",
    )
    ui.takeaway(
        "The only clue to the order of the reads is that neighbours overlap. Ordering reads by their "
        "overlaps is the travelling salesman problem solved in the rest of this application."
    )
    st.caption(
        "Continue with the page Learn (section Genome assembly) for the toy example and the objective."
    )
    ui.sources("sanger1977seq", "fleischmann1995", "lander1988")

else:
    st.subheader("Check yourself: mutations")
    ui.quiz("mutation")
