"""Page 1, Learn: the story in plain words, a toy example you can change, the objective, the glossary."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import content, figs  # noqa: E402

st.title("Learn: putting a genome back together")
st.caption("Everything on this page is computed live; no internet is needed.")

story, toy_tab, objective, glossary = st.tabs(
    ["The story", "Toy example", "The objective", "Glossary"]
)

with story:
    st.subheader("A shredded book")
    st.markdown(
        """
Imagine a book in which only four letters exist: **A, C, G and T**. That is a **genome**: the
instruction book of a living thing, for example a small virus. A sequencing machine cannot read the
whole book at once. It reads **many short random pieces**, called **reads**, and many pieces overlap
with their neighbours, like the torn edges of shredded pages.

The job: **put the pieces back in the right order**. Two pieces fit when the *end of one equals the
beginning of the other*. The longer that matching edge (the **overlap**), the more letters we save by
gluing the pieces, and the more likely they really were neighbours.
"""
    )
    st.subheader("From reads to a travelling salesman")
    st.markdown(
        """
1. **Reads** come from the genome with a few typing mistakes (read errors).
2. For every ordered pair of reads we compute the **overlap**. This gives a table that is *not*
   symmetric: A can fit before B without B fitting before A.
3. A **solution** is an *order* of all reads. Its **fitness F** is the total number of letters saved
   by gluing neighbours in that order. **Bigger is better.**
4. Finding the best order is the **travelling salesman problem** (visit every read once, travel as
   little as possible), in its asymmetric form (**ATSP**). The proof of this encoding is in
   `docs/OBJECTIVE_FUNCTION.md`.
5. Because the number of orders explodes (see below), we use clever search methods instead of trying
   them all: **greedy**, **tabu search**, **LKH-3** and **ALNS**, where a *bandit* learns which of
   several repair tools to use.
"""
    )
    st.subheader("Why we cannot try every order")
    rows = []
    for n in (5, 10, 20, 50, 100, 500):
        digits = len(str(math.factorial(n)))
        rows.append(
            {
                "reads n": n,
                "number of orders n!": f"{math.factorial(n):,}"
                if n <= 20
                else f"a number with {digits} digits",
            }
        )
    ui.table(rows)
    years = math.factorial(20) / 1e9 / (365.25 * 24 * 3600)
    st.caption(
        f"A computer that checks a billion orders per second would need about {years:.0f} years for 20 reads."
    )

with toy_tab:
    st.subheader("Toy example: three reads")
    st.markdown("Change the reads (letters A, C, G, T; 2 to 7 reads, at most 40 letters each).")
    left, right = st.columns([3, 1])
    text = left.text_input("Reads, separated by spaces", " ".join(content.TOY_READS))
    tau = right.number_input("Minimum overlap tau", 1, 10, 1, 1)
    reads = content.clean_reads(text)
    try:
        toy = content.make_toy(reads, int(tau))
    except ValueError as error:
        st.error(str(error))
        st.stop()
    for note in toy.notes:
        st.info(note)
    inst = toy.instance
    rows = content.all_orders(toy)
    best = rows[0]
    order_text = st.text_input(
        "Your order (every letter once, e.g. " + "".join(toy.labels) + ")", "".join(toy.labels)
    )
    try:
        mine = content.order_from_letters(toy, order_text)
    except ValueError as error:
        st.warning(str(error))
        mine = None
    c1, c2 = st.columns(2)
    with c1:
        ui.show_figure(
            figs.toy_graph(inst, toy.labels, None if mine is None else mine.tolist()),
            "toy_overlap_graph",
            caption="Red arrows: the glue steps of your order.",
        )
    with c2:
        ui.show_figure(
            figs.matrix_figure(inst.ov, toy.labels, "Overlap matrix", "letters saved"),
            "toy_overlap_matrix",
        )
    if mine is not None:
        f, length = inst.fitness(mine), inst.merged_length(mine)
        merged = " + ".join(inst.assemble(mine).sequences)
        m1, m2, m3 = st.columns(3)
        m1.metric(
            "Fitness F of your order",
            f,
            delta=int(f) - int(best["F"]),
            help="Compared with the best order",
        )
        m2.metric("Length L of the merged text", length)
        m3.metric("Best possible F", int(best["F"]))
        st.code(merged, language=None)
        st.caption(
            "Pieces joined by + are separate contigs: where neighbours do not overlap the chain breaks."
        )
    st.markdown(
        f"**All {len(rows)} orders** (best first). The total length of the reads is {int(inst.lengths.sum())}, and L = that total - F."
    )
    ui.table(rows)
    with st.expander("Greedy assembly, one decision at a time"):
        steps, greedy_order = content.greedy_steps(inst.ov)
        st.markdown(
            "Greedy sorts all overlaps from large to small and accepts an overlap when the first read "
            "has no successor yet, the second read no predecessor yet, and no cycle is created."
        )
        ui.table(
            [
                {
                    "overlap": f"{toy.labels[s['from']]} -> {toy.labels[s['to']]}",
                    "letters saved": s["overlap"],
                    "decision": s["verdict"],
                }
                for s in steps
            ]
        )
        names = "".join(toy.labels[i] for i in greedy_order)
        st.success(
            f"Greedy result: {names}, F = {inst.fitness(greedy_order)} (best possible: {int(best['F'])})"
        )

with objective:
    st.subheader("The objective function")
    st.markdown(
        "Put the reads in an order, glue every read to the one before it, and add up the letters saved. "
        "That total is the **fitness F**."
    )
    st.latex(r"F(\pi)=\sum_{k=1}^{n-1}\operatorname{ov}(\pi_k,\pi_{k+1}) \qquad\text{(maximise)}")
    st.latex(
        r"L(\pi)=\Big(\sum_{i=1}^{n}\ell_i\Big)-F(\pi) \qquad\text{(the merged text is this long)}"
    )
    st.markdown(
        "`ov(i, j)` is the length of the longest overlap where the end of read *i* equals the beginning of "
        "read *j* (at least tau letters, else 0). With read errors we allow a few mismatches in an overlap."
    )
    st.markdown(
        "**The same problem as a travelling salesman.** Add a dummy node 0 and use the cost"
    )
    st.latex(r"c(i,j)=\ell_j-\operatorname{ov}(i,j),\quad c(0,j)=\ell_j,\quad c(i,0)=0")
    st.markdown(
        "A tour of the nodes costs exactly L, so *minimising the tour = minimising the merged length = "
        "maximising F*. The costs are not symmetric, so it is an **ATSP**."
    )
    st.info(
        "What F does not know: a high F does not guarantee the real genome. With read errors, the best "
        "order we find can have *higher* F than the true order (decision D-036), which is why we also "
        "measure identity, genome fraction and contigs against the real genome."
    )

with glossary:
    st.subheader("Glossary")
    st.caption("Read from docs/glossary.yaml, the same source as the report and the slides.")
    terms = content.load_glossary()
    query = st.text_input("Search the glossary", placeholder="e.g. coverage, bandit, tabu")
    hits = content.search_glossary(terms, query)
    st.write(f"{len(hits)} of {len(terms)} terms")
    for term in hits:
        with st.expander(f"{term['term']}  (Phase {term.get('phase', '?')})"):
            st.markdown(f"**In plain words.** {str(term['plain']).strip()}")
            st.markdown(f"**Precisely.** {str(term['precise']).strip()}")
