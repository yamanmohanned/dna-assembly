"""Presentation: the whole lecture as a full-screen slide deck with live 3D models and animations."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import lecture, studio  # noqa: E402


@st.cache_data(show_spinner="Building the lecture from the project data...", max_entries=4)
def _deck(author: str, embedded: bool, stamp: str) -> str:
    """The deck as HTML; ``stamp`` makes the cache notice an edited script or style sheet."""
    return lecture.html(author, embedded=embedded)


st.title("Presentation", icon=":material/slideshow:")
left, right = st.columns([3, 2], vertical_alignment="bottom")
left.markdown(
    "The lecture **from the double helix to the assembled genome**, in six parts. Click the deck once, "
    "then present with the keyboard."
)
author = right.text_input(
    "Name on the title slide (optional)", key="lecture-author", max_chars=80,
    placeholder="Your name, course, date",
)  # fmt: skip
st.markdown(
    ":blue-badge[→ or Space: next] :blue-badge[←: back] :blue-badge[F: full screen] "
    ":blue-badge[O: overview] :blue-badge[P: play or pause an animation]"
)
stamp = studio.assets_stamp()
ui.studio(_deck(author.strip(), True, stamp), 760, "Lecture slide deck")
st.caption(
    "On a slide with an animation, → plays one step at a time before moving on. 3D models can be turned "
    "with the mouse during the talk. Every number on the slides is computed from the project data or "
    "quoted with its source; the last slide lists the references."
)
c1, c2 = st.columns([1, 3], vertical_alignment="center")
c1.download_button(
    "Lecture as one HTML file",
    data=lambda: _deck(author.strip(), False, stamp).encode("utf-8"),
    file_name="dna_lecture.html",
    mime="text/html",
    icon=":material/download:",
    on_click="ignore",
)
c2.caption(
    "The file needs neither this application nor the internet: open it by double-click in any modern "
    "browser and press F for full screen. Keep it as the fallback for the day of the presentation."
)
with st.expander("Outline of the lecture", icon=":material/list:"):
    deck = lecture.build()
    ui.table(
        [
            {
                "Slide": i + 1,
                "Part": slide["section"],
                "Title": str(slide["title"]).replace("<br>", " "),
            }
            for i, slide in enumerate(deck["slides"])
        ]
    )
