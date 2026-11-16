"""DNA structure: the double helix in 3D (real atoms), the nucleotide, base pairs, measured numbers."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import biology, studio  # noqa: E402
from genomeasm.appkit.paths import GENOMES  # noqa: E402
from genomeasm.molbio import diagrams, library, sequence  # noqa: E402

st.title("DNA structure", icon=":material/genetics:")
st.markdown(
    "**Deoxyribonucleic acid (DNA)** is a long chain of four kinds of nucleotides. The **order** of "
    "the nucleotides is the hereditary information of a cell. Two such chains wind around each "
    "other as the **double helix**."
)
measured = library.measured_helix()
facts = library.facts()
section = ui.section_picker(
    "Section",
    [
        "Double helix",
        "Nucleotide",
        "Base pairs",
        "Helix in numbers",
        "Chargaff's rule",
        "Check yourself",
    ],
    "structure-section",
)

if section == "Double helix":
    st.subheader("Two antiparallel strands, wound into a right-handed helix")
    ui.key_facts(
        [
            (
                "Rise per base pair",
                f"{measured['rise_mean']:.2f} Å",
                "Measured by this program on PDB 1BNA (mean of 11 steps)",
            ),
            (
                "Base pairs per turn",
                f"{measured['bp_per_turn']:.1f}",
                "Measured on PDB 1BNA; in solution 10.4 ± 0.1 (Wang 1979)",
            ),
            (
                "Hydrogen bonds",
                "2 in A·T, 3 in G·C",
                f"All {int(measured['n_hbonds'])} bonds of the {int(measured['n_pairs'])} pairs of 1BNA were found from the coordinates",
            ),
            ui.fact_item("diameter"),
        ]
    )
    model_kind = st.segmented_control(
        "Model",
        ["Experimental structure (real atoms)", "Teaching model of any sequence"],
        default="Experimental structure (real atoms)",
        key="structure-model",
        required=True,
    )
    if model_kind == "Experimental structure (real atoms)":
        structure = library.structure_1bna()
        ui.studio(
            studio.viewer_html(
                biology.scenes()["atomic"], height=600, spin=True, colour="strand",
                note="PDB 1BNA, X-ray diffraction, 1.9 Å; hydrogen atoms are not resolved.",
            ),
            690,
            "Interactive 3D model of the DNA dodecamer 1BNA",
        )  # fmt: skip
        st.caption(
            f"**{structure.title.capitalize()}** PDB entry {structure.pdb_id}: {len(structure.atoms)} DNA atoms in two "
            f"chains of sequence {structure.sequence('A')}, {structure.experiment.lower()} at {structure.resolution} Å "
            "(Drew et al. 1981). The file is in `data/raw` with its SHA-256 checksum."
        )
    else:
        c1, c2 = st.columns([3, 1])
        accession = c1.selectbox(
            "Genome",
            list(GENOMES),
            format_func=lambda a: f"{a}: {GENOMES[a]}",
            key="structure-genome",
        )
        last = len(biology.genome(accession)) - biology.SCHEMATIC_LETTERS + 1
        first = int(c2.number_input("First letter", 1, last, 1, 1, key="structure-first"))
        scene = biology.schematic(accession, first)
        ui.studio(
            studio.viewer_html(
                scene, height=600, spin=True, colour="base",
                note="Teaching model: one bead per phosphate, sugar and base.",
            ),
            690,
            "Interactive 3D teaching model of a double helix",
        )  # fmt: skip
        piece = biology.genome(accession)[first - 1 : first - 1 + biology.SCHEMATIC_LETTERS]
        st.code(
            f"5′ {piece} 3′   strand 1 (letters {first:,} to {first + len(piece) - 1:,})\n"
            f"3′ {sequence.complement(piece)} 5′   strand 2",
            language=None,
        )
        st.caption(
            f"**A simplified model, not a measurement.** The beads are placed with the average backbone geometry "
            f"measured on 1BNA and stacked with the published rise ({scene['rise']} Å) and helical repeat in solution "
            f"({scene['bp_per_turn']} base pairs per turn, Wang 1979). This stretch is held together by "
            f"{sequence.hydrogen_bonds(piece)} hydrogen bonds."
        )
    with st.expander("What to notice in the model", icon=":material/visibility:"):
        st.markdown(
            "- **Antiparallel strands.** Find the labels 5′ and 3′: at each end of the helix one strand ends "
            "with 5′ and the other with 3′.\n"
            "- **Backbone outside, bases inside.** Colour by *Part* to see the phosphates and sugars on the "
            "outside and the flat bases stacked in the middle.\n"
            "- **Two grooves.** The backbones of a base pair are not opposite each other, so the helix has a "
            "wide **major groove** and a narrow **minor groove**. Proteins that recognise a sequence read the "
            "edges of the bases in the grooves.\n"
            "- **Right-handed.** Follow one backbone upwards: it turns like a normal screw.\n"
            "- **Click a nucleotide** to see its partner and the hydrogen bonds between them."
        )
    ui.takeaway(
        "The two strands are complementary and antiparallel: each one fully determines the other. "
        "That is why DNA can be copied, repaired and read."
    )
    ui.sources("watson1953", "franklin1953", "wilkins1953", "drew1981", "wang1979")

elif section == "Nucleotide":
    st.subheader("The building block: phosphate, sugar, base")
    left, right = st.columns(2)
    with left:
        ui.diagram(
            diagrams.nucleotide_svg(),
            "Schematic drawing of a nucleotide",
            "Schematic. The primes (1′ to 5′) number the carbons of the sugar.",
        )
        st.markdown(
            "- The **base** is attached to carbon **1′** of the sugar.\n"
            "- The **phosphate** is attached to carbon **5′**.\n"
            "- Carbon **3′** carries a free **OH** group: the next nucleotide is joined there by a "
            "**phosphodiester bond**.\n"
            "- A strand therefore has two different ends and a **direction, written 5′ → 3′**. "
            "DNA polymerase and RNA polymerase can only add to the 3′ end."
        )
    with right:
        ui.studio(
            studio.viewer_html(
                biology.scenes()["nucleotide"], height=430, colour="component", atomLabels=True,
                selectable=False, hint="Drag to rotate · a real nucleotide of PDB 1BNA",
            ),
            520,
            "3D model of one deoxyadenosine nucleotide",
        )  # fmt: skip
        st.caption(
            "A real deoxyadenosine nucleotide (residue 5 of chain A of 1BNA) with the atom names of the PDB "
            "file. Hydrogen atoms are not resolved at 1.9 Å."
        )
    ui.takeaway(
        "The sugar-phosphate backbone is identical along the whole molecule. Only the bases differ, and "
        "their order is the information."
    )
    ui.sources("drew1981", "pdb1bna", "alberts2022")

elif section == "Base pairs":
    st.subheader("A pairs with T, G pairs with C")
    st.markdown(
        "A **purine** (adenine or guanine, two fused rings) always faces a **pyrimidine** (thymine or "
        "cytosine, one ring). They are held together by **hydrogen bonds** between the atoms on their "
        "facing edges."
    )
    ui.key_facts(
        [
            ("A·T pair", "2 hydrogen bonds", "N6-H···O4 and N3-H···N1"),
            ("G·C pair", "3 hydrogen bonds", "N4-H···O6, N1-H···N3 and N2-H···O2"),
            (
                "Bond length",
                f"{measured['hbond_mean']:.2f} Å on average",
                f"Measured on PDB 1BNA: {measured['hbond_min']:.2f} to {measured['hbond_max']:.2f} Å between the two heavy atoms",
            ),
            (
                "Width of a pair",
                f"{measured['c1_distance_mean']:.1f} Å",
                "Distance C1′ to C1′, the same for A·T and G·C (measured on PDB 1BNA)",
            ),
        ]
    )
    left, right = st.columns(2)
    for column, name, title in (
        (left, "pair_at", "Adenine · thymine"),
        (right, "pair_gc", "Guanine · cytosine"),
    ):
        with column:
            st.markdown(f"**{title}**")
            ui.studio(
                studio.viewer_html(
                    biology.scenes()[name], height=380, colour="element", distanceLabels=True,
                    atomLabels=True, selectable=False, toolbar=False,
                    hint="Drag to rotate · distances in ångström",
                ),
                430,
                f"3D model of a real {title} base pair",
            )  # fmt: skip
    st.caption(
        "Two real base pairs of PDB 1BNA seen along the helix axis, with sugar and phosphate. Dotted lines: "
        "hydrogen bonds; the numbers are distances measured between the bonded atoms."
    )
    ui.takeaway(
        "Because both pairs have the same width, the helix is equally thick everywhere, whatever its "
        "sequence. G·C pairs, with three bonds, hold the strands together more strongly than A·T pairs."
    )
    ui.sources("watson1953", "pdb1bna")

elif section == "Helix in numbers":
    st.subheader("Measured here, compared with the literature")
    st.markdown(
        "The first column of numbers is **measured by this program** on the atom coordinates of PDB 1BNA; "
        "the second is the **published** value with its source."
    )
    rows = biology.helix_rows()
    ui.table(
        [
            {
                "Quantity": r["quantity"],
                "Measured on PDB 1BNA": r["measured"],
                "Published": r["published"],
            }
            for r in rows
        ]
    )
    with st.expander("How the helix is measured", icon=":material/straighten:"):
        st.markdown(
            "1. **Base pairs** are found from the coordinates: a purine N1 atom within 3.5 Å of a pyrimidine "
            "N3 atom of the other chain.\n"
            "2. The **centre** of a pair is the midpoint of the purine C8 and the pyrimidine C6 atom, which "
            "lies close to the helix axis in B-DNA.\n"
            "3. The **axis** is the straight line that fits the twelve centres best.\n"
            "4. **Rise** is the distance between neighbouring centres along the axis; **twist** is the angle "
            "between the C8/C6 lines of neighbouring pairs seen along the axis. A positive twist means a "
            "right-handed helix.\n\n"
            "A crystal structure is slightly bent (Drew et al. report 19 degrees for 1BNA), so a straight "
            "axis is an approximation; this is why standard deviations are given. The code is in "
            "`genomeasm/molbio/structure.py`."
        )
    ui.takeaway(
        "The measurements on real coordinates agree with the published structure of B-DNA: about 3.4 Å "
        "per base pair and about ten base pairs per right-handed turn."
    )
    ui.sources("drew1981", "watson1953", "wang1979", "pdb1bna")

elif section == "Chargaff's rule":
    st.subheader("As much A as T, as much G as C, if there are two strands")
    st.markdown(
        "If A always pairs with T and G with C, a **double-stranded** DNA must contain equal amounts of A "
        "and T and equal amounts of G and C. Chargaff found this regularity in 1950, before the double "
        "helix was known. The two genomes of this project put it to the test."
    )
    rows = biology.composition_rows()
    ui.diagram(
        diagrams.composition_svg(rows),
        "Base composition of phiX174 and lambda",
        "Counted from the sequences in data/raw. The axis starts at zero.",
    )
    table = []
    for name, counts in rows:
        total = sum(counts.values())
        table.append(
            {
                "Genome": name,
                "A %": round(100 * counts["A"] / total, 2),
                "T %": round(100 * counts["T"] / total, 2),
                "G %": round(100 * counts["G"] / total, 2),
                "C %": round(100 * counts["C"] / total, 2),
                "A / T": round(counts["A"] / counts["T"], 3),
                "G / C": round(counts["G"] / counts["C"], 3),
            }
        )
    ui.table(table)
    st.markdown(
        "- **phiX174** breaks the rule (A / T far from 1), and correctly so: inside the virus particle its "
        "genome is a **single strand**, so there is no partner strand that would balance A against T.\n"
        "- **lambda** is double-stranded. The file holds one of its two strands; counted over **both** "
        "strands A equals T exactly, by construction. That the ratios are already close to 1 on a single "
        "strand is an observed regularity of long genomes, not a consequence of base pairing."
    )
    ui.takeaway("Unequal amounts of A and T are the signature of a single-stranded genome.")
    ui.sources("chargaff1950", "refseq_phix174", "refseq_lambda")

else:
    st.subheader("Check yourself: structure")
    ui.quiz("structure")
