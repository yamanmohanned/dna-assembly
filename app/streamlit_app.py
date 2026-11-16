"""DNA, from the double helix to genome assembly as a travelling salesman problem: the application.

Start it with ``genomeasm app run`` (or ``streamlit run app/streamlit_app.py``). It works offline:
the genomes come from ``data/raw``, the precomputed results from ``results/phase5``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

st.set_page_config(
    page_title="DNA: from the double helix to genome assembly",
    page_icon=":material/genetics:",
    layout="wide",
    initial_sidebar_state="expanded",
)

pages = {
    "Start": [
        st.Page("views/guide.py", title="Guide", icon=":material/menu_book:", default=True),
        st.Page("views/presentation.py", title="Presentation", icon=":material/slideshow:"),
    ],
    "Molecular biology": [
        st.Page("views/dna_structure.py", title="DNA structure", icon=":material/genetics:"),
        st.Page(
            "views/central_dogma.py", title="DNA to protein", icon=":material/conversion_path:"
        ),
        st.Page("views/genome.py", title="Genes and mutations", icon=":material/account_tree:"),
    ],
    "Genome assembly": [
        st.Page("views/learn.py", title="Learn", icon=":material/school:"),
        st.Page("views/data.py", title="Data", icon=":material/database:"),
        st.Page("views/simulator.py", title="Simulator", icon=":material/tune:"),
        st.Page("views/overlap_graph.py", title="Overlap Graph", icon=":material/hub:"),
        st.Page("views/solver_lab.py", title="Solver Lab", icon=":material/science:"),
        st.Page("views/experiments.py", title="Experiments", icon=":material/bar_chart:"),
        st.Page(
            "views/assembly_viewer.py", title="Assembly Viewer", icon=":material/scatter_plot:"
        ),
    ],
    "Sources": [
        st.Page("views/references.py", title="References", icon=":material/library_books:"),
        st.Page("views/reproducibility.py", title="Reproducibility", icon=":material/verified:"),
    ],
}
ui.force_left_to_right()
ui.logo()
navigation = st.navigation(pages, expanded=True)
st.sidebar.caption(
    "Offline application. Genomes: NCBI RefSeq (phiX174, lambda); structure: PDB 1BNA. Live runs "
    "end in under 30 s; larger runs are shown from precomputed results and labelled."
)
navigation.run()
