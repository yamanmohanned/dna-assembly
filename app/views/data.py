"""Page 2, Data: the genome cards, GC content and the data-loading pipeline."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import figs, lab, repro  # noqa: E402
from genomeasm.appkit.paths import GENOMES, PROCESSED, RAW  # noqa: E402

st.title("Data: two real genomes, checked and stored offline")
st.markdown(
    "All sequences are **real NCBI sequences**. They were downloaded once, checked, fingerprinted and "
    "stored in `data/raw/`; this application only reads the stored copy, so it works without internet."
)

accession = st.radio(
    "Genome", list(GENOMES), format_func=lambda a: f"{a}: {GENOMES[a]}", horizontal=True
)
card = json.loads((PROCESSED / f"{accession}.card.json").read_text(encoding="utf-8"))

c1, c2, c3, c4 = st.columns(4)
c1.metric("Length (letters)", f"{card['length_bp']:,}")
c2.metric("GC content", f"{card['gc_percent']:.1f}%", help="Share of the letters G and C")
c3.metric("Molecule", card["molecule_type"].upper() + f" ({card['strand_reported_by_ncbi']})")
c4.metric("Topology (NCBI)", card["topology_reported_by_ncbi"])
st.markdown(
    f"**{card['organism']}**, version {card['accession_version']}, NCBI {card['ncbi_source_database']} "
    f"record, last updated {card['ncbi_update_date']}. We treat every genome as one *linear* string "
    "(decision D-017), whatever its real topology."
)

left, right = st.columns([2, 1])
with left:
    window = st.slider("GC window (letters)", 50, 1000, 200, 50)
    ui.show_figure(
        figs.gc_profile(lab.genome_sequence(accession), window),
        f"gc_profile_{accession}",
        caption=f"Source: {card['file']} (sha256 {card['file_sha256'][:12]}...)",
    )
with right:
    ui.show_figure(figs.base_counts(card), f"base_counts_{accession}")

with st.expander("The first lines of the FASTA file"):
    lines = (RAW / card["file"]).read_text(encoding="ascii").splitlines()[:4]
    st.code("\n".join(lines), language=None)
    st.caption("FASTA: one header line starting with >, then the letters in lines of 70.")

st.subheader("The data-loading pipeline")
ui.show_figure(
    figs.pipeline_diagram(),
    "data_loading_pipeline",
    caption=(
        "Orange: runs once with internet (genomeasm data fetch). Green: runs on every start, offline. "
        "Details: docs/DATA_LOADING.md."
    ),
)
st.markdown(
    "We copy the genome **once** from the NCBI library, check the copy against the library's own "
    "facts (name, kind of molecule, length) and save it with a fingerprint (SHA-256). From then on the "
    "program reads only the saved copy and checks the fingerprint first, so a single changed letter is "
    "noticed."
)
if st.button("Check the fingerprints of my stored genomes now", icon=":material/verified:"):
    rows = repro.verify_checksums()
    ui.table(
        [
            {
                "file": r["file"],
                "fingerprint matches": "yes" if r["ok"] else "NO",
                "sha256": r["actual"],
            }
            for r in rows
        ]
    )
    (st.success if all(r["ok"] for r in rows) else st.error)(
        "All stored genomes are exactly as downloaded."
        if all(r["ok"] for r in rows)
        else "A stored genome differs from its fingerprint."
    )
