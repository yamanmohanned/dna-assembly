"""Page 0, Guide: everything a newcomer needs, the presentation script and the glossary."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import content  # noqa: E402
from genomeasm.appkit.paths import DOCS  # noqa: E402

st.title("Guide: start here")
ui.banner()
st.info(
    "No background needed. Read sections 1 to 6 to understand the project, section 7 to use the "
    "application, section 8 to present it, section 9 for questions. The glossary explains every term.",
    icon=":material/menu_book:",
)
sections = content.guide_sections()
tabs = st.tabs([title for title, _ in sections] + ["Glossary"])
for tab, (_, body) in zip(tabs, sections, strict=False):
    with tab:
        st.markdown(body)
with tabs[-1]:
    terms = content.load_glossary()
    query = st.text_input(
        "Search the glossary", key="guide-glossary", placeholder="e.g. overlap, bandit"
    )
    hits = content.search_glossary(terms, query)
    st.write(f"{len(hits)} of {len(terms)} terms")
    for term in hits:
        with st.expander(str(term["term"])):
            st.markdown(f"**In plain words.** {str(term['plain']).strip()}")
            st.markdown(f"**Precisely.** {str(term['precise']).strip()}")
c1, c2, _ = st.columns([1, 1, 4])
c1.download_button(
    "Guide (Markdown)", (DOCS / "GUIDE.md").read_bytes(), "GUIDE.md", "text/markdown",
    icon=":material/download:", on_click="ignore",
)  # fmt: skip
c2.download_button(
    "Glossary (Markdown)", content.glossary_markdown(terms).encode("utf-8"), "GLOSSARY.md",
    "text/markdown", icon=":material/download:", on_click="ignore",
)  # fmt: skip
