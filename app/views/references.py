"""References: every source of the molecular-biology pages and the lecture, with how it was checked."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.molbio import library  # noqa: E402

st.title("References", icon=":material/library_books:")
refs = library.references()
facts = library.facts()
checked = sum(1 for r in refs.values() if r["check"] == "crossref")
quoted = sum(1 for f in facts.values() if f["check"] == "abstract")
open_facts = sum(1 for f in facts.values() if f["check"] == "textbook")
ui.key_facts(
    [
        ("References", str(len(refs)), "Journal articles, one textbook and the data sources"),
        (
            "Compared with Crossref",
            str(checked),
            "Title, journal, volume, pages, year and authors of every journal article",
        ),
        (
            "Numbers read at the source",
            str(quoted),
            "Read in the abstract of the paper (PubMed record)",
        ),
        (
            "Still to verify",
            str(open_facts),
            "Textbook values and paper conclusions that were not re-read at the source",
        ),
    ]
)
st.markdown(
    "**What is established and what is simplified.** Statements on the pages are either *measured* by this "
    "program from the data in `data/raw`, *quoted* from a source listed here, or part of a *schematic* "
    "drawing, which is labelled as such. Nothing is cited from memory without saying so."
)
view = ui.section_picker(
    "View", ["References", "Quoted facts", "How things were checked"], "refs-view"
)

if view == "References":
    query = st.text_input(
        "Search the references", placeholder="e.g. helix, Sanger, 1977", key="refs-search"
    )
    words = [w for w in query.lower().split() if w]
    shown = 0
    for topic in library.TOPIC_ORDER:
        group = [
            r for r in refs.values()
            if r["topic"] == topic and all(w in library.reference_line(r).lower() for w in words)
        ]  # fmt: skip
        if not group:
            continue
        st.subheader(topic.capitalize())
        for ref in sorted(group, key=lambda r: (str(r["authors"][0]), int(r["year"] or 0))):
            shown += 1
            colour = {"crossref": "green", "file": "blue", "software": "blue"}.get(
                str(ref["check"]), "orange"
            )
            line = library.reference_line(ref)
            if ref.get("doi"):
                line = line.replace(
                    f"doi:{ref['doi']}", f"[doi:{ref['doi']}](https://doi.org/{ref['doi']})"
                )
            st.markdown(f"{line}  \n:{colour}-badge[{library.CHECK_LABELS[str(ref['check'])]}]")
            if ref.get("note"):
                st.caption(str(ref["note"]))
    st.caption(
        f"{shown} of {len(refs)} references shown. The DOI links need the internet; everything else works offline."
    )
    st.download_button(
        "References (text)",
        data="\n".join(library.reference_line(r) for r in refs.values()).encode("utf-8"),
        file_name="references.txt",
        mime="text/plain",
        icon=":material/download:",
        on_click="ignore",
    )

elif view == "Quoted facts":
    st.markdown(
        "Every number or statement taken from the literature, with its source and the way it was checked. "
        "Values measured by the program itself are not in this table; they are shown where they are computed."
    )
    ui.table(
        [
            {
                "Topic": f["topic"],
                "Fact": f["label"],
                "Value": f["value"],
                "Source": library.cite(str(f["source"])),
                "Checked": library.CHECK_LABELS[str(f["check"])],
                "Statement": f["statement"],
            }
            for f in facts.values()
        ],
        height=560,
    )

else:
    st.markdown(
        "| Kind of content | How it was checked |\n|---|---|\n"
        "| Bibliographic data of journal articles | Compared with the Crossref record of each DOI on 2026-11-14 |\n"
        "| Numbers quoted from papers | Read in the abstract of the PubMed record on 2026-11-14, where the abstract states them |\n"
        "| Textbook values (fork speed, fidelity, Okazaki fragment length, 46 chromosomes) | **Not** re-read at the source: marked TODO-VERIFY |\n"
        "| Helix geometry, hydrogen bonds | Measured by this program on PDB 1BNA; tests compare them with the published ranges |\n"
        "| Genetic code | Built from NCBI translation table 1; a test compares it with Biopython |\n"
        "| Genes of phiX174 | Coordinates read from the NCBI feature table; a test checks each gene against the sequence |\n"
        "| Genome sequences and the structure file | SHA-256 checksums in `data/raw/SHA256SUMS`, shown on the page Reproducibility |\n"
        "| Results of the experiments | Read from `results/` by the same loader as the report |"
    )
    m = library.measured_helix()
    st.caption(
        f"Example of a measured value: rise {m['rise_mean']:.2f} ± {m['rise_sd']:.2f} Å per base pair and "
        f"{m['bp_per_turn']:.1f} base pairs per turn on PDB 1BNA; Drew et al. (1981) report 9.8 base pairs "
        "per turn for the centre of the same crystal."
    )
