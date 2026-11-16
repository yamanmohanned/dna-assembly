"""The lecture: one slide deck from the double helix to the assembled genome.

``build()`` returns the deck as data (slides, 3D scenes, animation scenarios, references) and
``html()`` the self-contained document shown by the Presentation page and written by
``genomeasm app lecture``. Every number on a slide is computed here from the project's data
(``results/``, the genomes, PDB 1BNA) or quoted from ``docs/molbio_facts.yaml`` with its source;
slide text never repeats a number by hand.
"""

from __future__ import annotations

import base64
from html import escape
from pathlib import Path
from typing import Any

import numpy as np

from genomeasm.appkit import biology, content, figs, studio
from genomeasm.appkit.paths import FIGURES
from genomeasm.molbio import diagrams, library, sequence

TITLE = "DNA: from the double helix to the assembled genome"
SECTIONS = (
    ("Structure", "What DNA is made of and why it is a double helix"),
    ("Genomes", "Genes, chromosomes and a whole genome on one slide"),
    ("Replication", "How one helix becomes two identical ones"),
    ("Expression", "From DNA to RNA to protein: the genetic code"),
    ("Mutation and biotechnology", "When the text changes, and how we read and edit it"),
    ("Sequencing and assembly", "This project: putting a genome back together"),
)


def _cites(*keys: str) -> list[str]:
    return list(dict.fromkeys(library.cite(k) for k in keys))  # in order, without repeats


def _fact(key: str, tone: str = "") -> dict[str, str]:
    fact = library.facts()[key]
    source = library.cite(str(fact["source"]))
    if fact["check"] == "textbook":
        source += " (textbook value)"
    return {
        "value": str(fact["value"]),
        "label": str(fact["label"]),
        "source": source,
        "tone": tone,
    }


def _data_uri(path: Path, mime: str) -> str:
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def _svg_uri(svg: bytes) -> str:
    return "data:image/svg+xml;base64," + base64.b64encode(svg).decode("ascii")


def _table(headers: list[str], rows: list[list[str]], numeric: tuple[int, ...] = ()) -> str:
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join(
        "<tr>"
        + "".join(
            f'<td class="num">{c}</td>' if i in numeric else f"<td>{c}</td>"
            for i, c in enumerate(row)
        )  # fmt: skip
        + "</tr>"
        for row in rows
    )
    return f'<table class="deck-table"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def _results() -> dict[str, Any]:
    """Numbers of the experiments, read from ``results/`` (the same loader as the report)."""
    from genomeasm.docgen.facts import Facts, name, p_text

    f = Facts.load()
    table = f.by_algorithm
    order = ["greedy", "tabu", "alns_random", "alns_ucb1", "alns_ucb1_rt", "lkh3"]
    rows = [
        [name(a), f"{table.loc[a, 'mean_gap']:.2f}%", f"{table.loc[a, 'bio_n_contigs']:.0f}",
         f"{100 * table.loc[a, 'bio_identity']:.2f}%"]
        for a in order
    ]  # fmt: skip
    ucb, eps = f.verdict_ucb, f.verdict_eps
    post = f.p5.info["post_hoc_ucb1_tie_breaking"][0]
    return {
        "n_sets": f.n_sets, "n_settings": f.n_settings, "seeds": f.seeds, "budget": f.budget,
        "rows": rows, "ucb_p": p_text(float(ucb["p_holm"])), "ucb_a12": float(ucb["a12"]),
        "eps_a12": float(eps["a12"]), "post_a12": float(post["a12"]),
        "post_p": p_text(float(post["p_holm"])), "lkh_never_beaten": bool(f.lkh_never_beaten),
    }  # fmt: skip


def _toy_figure() -> tuple[str, dict[str, Any]]:
    toy = content.make_toy(list(content.TOY_READS), 1)
    best = content.all_orders(toy)[0]
    order = content.order_from_letters(toy, str(best["order"]))
    figure = figs.toy_graph(toy.instance, tuple(toy.labels), tuple(int(i) for i in order))
    info = {
        "reads": list(content.TOY_READS), "order": best["order"], "f": int(best["F"]),
        "l": int(best["L"]), "merged": best["merged"],
        "total": int(np.sum(toy.instance.lengths)),
    }  # fmt: skip
    return _svg_uri(figs.figure_bytes(figure, "svg")), info


def _references_html() -> list[str]:
    """The reference list as two slides of equal length (by topic; a topic may run over)."""
    refs = library.references()
    entries: list[tuple[str, str]] = []  # (topic, one reference as HTML)
    for topic in library.TOPIC_ORDER:
        group = [r for r in refs.values() if r["topic"] == topic]
        for ref in sorted(group, key=lambda r: (str(r["authors"][0]), int(r["year"] or 0))):
            mark = (' <span class="unchecked">[not checked: TODO-VERIFY]</span>'
                    if ref["check"] == "not-checked" else "")  # fmt: skip
            entries.append((topic, f"<p>{escape(library.reference_line(ref))}{mark}</p>"))
    half = sum(len(text) + 60 for _, text in entries) / 2  # text length decides the height
    pages: list[list[str]] = [[], []]
    used, last = 0.0, ["", ""]
    for topic, text in entries:
        page = 0 if used + (len(text) + 60) / 2 <= half else 1
        if last[page] != topic:
            again = " (continued)" if topic in last else ""
            pages[page].append(f'<p class="ref-topic">{escape(topic)}{again}</p>')
            last[page] = topic
        pages[page].append(text)
        used += len(text) + 60
    return ["".join(page) for page in pages]


def build(author: str = "") -> dict[str, Any]:
    """The lecture as data. ``author`` is printed on the title slide when given."""
    facts = library.facts()
    m = library.measured_helix()
    phix, lam = biology.genome(biology.PHIX), biology.genome(biology.LAMBDA)
    genes = library.phix_genes()
    overlap = biology.overlap_summary()
    shotgun = biology.shotgun_example()
    gel_letters, gel_first = biology.sanger_window()
    res = _results()
    toy_uri, toy = _toy_figure()
    comp = dict(biology.composition_rows())
    px = sequence.base_composition(phix)
    lm = sequence.base_composition(lam)
    gene_j = library.gene("J")
    protein_j = sequence.translate(gene_j.sequence(phix))
    n_stop = len(sequence.STOP_CODONS)
    slides: list[dict[str, Any]] = []

    def add(section: str, layout: str, title: str, **fields: Any) -> None:
        number = next((i + 1 for i, s in enumerate(SECTIONS) if s[0] == section), 0)
        kicker = f"{number} · {section}" if number else section
        slides.append({"section": section, "layout": layout, "title": title, "kicker": kicker,
                       **fields})  # fmt: skip

    meta = "DNA fragment assembly as a travelling salesman problem<br>Genomes: NCBI RefSeq · "
    meta += "Structure: PDB 1BNA · every number is traceable to its source"
    if author:
        meta = escape(author) + "<br>" + meta
    add("Start", "title", "DNA: from the double helix<br>to the assembled genome",
        kicker="Molecular biology and computation",
        lead="Structure, replication and expression of DNA, and how a sequenced genome is put "
        "back together by solving a travelling salesman problem.",
        meta=meta,
        figure={"type": "viewer", "scene": "schematic",
                "options": {"spin": True, "labels": False, "selectable": False,
                            "colour": "base", "zoom": 1.12}})  # fmt: skip
    add("Start", "outline", "The lecture in six parts", kicker="Outline", items=[])

    # ---- 1 structure -----------------------------------------------------------------------------
    add("Structure", "split", "DNA stores hereditary information as a sequence",
        lead="<em class='term'>Deoxyribonucleic acid (DNA)</em> is a long chain of four kinds of "
        "<b>nucleotides</b>. The <b>order</b> of the nucleotides is the information a cell "
        "inherits and reads.",
        points=[
            "Four letters: <b>A</b>denine, <b>T</b>hymine, <b>G</b>uanine, <b>C</b>ytosine.",
            "Experiments of 1944 and 1952 showed that DNA, not protein, is the genetic material.",
            "In 1953 Watson and Crick proposed the <b>double helix</b>, built on the X-ray "
            "diffraction patterns of Franklin, Gosling and Wilkins.",
        ],
        facts=[{"value": "4", "label": "kinds of nucleotide", "source": "", "tone": ""},
               {"value": "2", "label": "strands, running in opposite directions", "source": "",
                "tone": "good"}],
        figure={"type": "viewer", "scene": "schematic",
                "options": {"spin": True, "colour": "base", "zoom": 1.05},
                "caption": f"Teaching model of letters 1 to {biology.SCHEMATIC_LETTERS} of the "
                "phiX174 genome: one bead per phosphate, sugar and base."},
        sources=_cites("avery1944", "hershey1952", "watson1953", "franklin1953", "wilkins1953"))  # fmt: skip
    add("Structure", "split", "A nucleotide has three parts",
        lead="Each nucleotide is a <b>phosphate</b> group, the sugar <b>deoxyribose</b> and one "
        "of four <b>bases</b>.",
        points=[
            "The carbons of the sugar are numbered 1′ to 5′. The base sits on <b>1′</b>, the "
            "phosphate on <b>5′</b>, a free OH group on <b>3′</b>.",
            "Nucleotides are joined 3′ to 5′ by <b>phosphodiester bonds</b>: this is the "
            "<b>sugar-phosphate backbone</b>.",
            "So every strand has a <b>direction</b>, written <b>5′ → 3′</b>. All enzymes that "
            "build DNA or RNA add to the 3′ end.",
        ],
        takeaway="The backbone is the same everywhere; <b>only the bases differ</b>, and they "
        "carry the information.",
        figure={"type": "html", "html": diagrams.nucleotide_svg(),
                "caption": "Schematic drawing; the next slide shows the real atoms."},
        cls="")  # fmt: skip
    n_atoms = len(biology.scenes()["nucleotide"]["atoms"]["el"])
    add("Structure", "split", "The same nucleotide, as measured by X-ray crystallography",
        lead="One real <b>deoxyadenosine nucleotide</b> from the crystal structure of a DNA "
        "double helix (PDB entry 1BNA, resolution 1.9 Å).",
        points=[
            f"{n_atoms} atoms other than hydrogen; hydrogen is too light to be seen at this "
            "resolution.",
            "<b>Orange</b>: phosphate group. <b>Violet</b>: the five-membered sugar ring. "
            "<b>Blue</b>: adenine, a flat system of two fused rings.",
            "The flat base is what allows bases to <b>stack</b> on each other inside the helix.",
        ],
        figure={"type": "viewer", "scene": "nucleotide",
                "options": {"colour": "component", "atomLabels": True, "selectable": False,
                            "style": "ballstick", "zoom": 1.0},
                "caption": "Atom names as in the PDB file. Drag to rotate."},
        sources=_cites("drew1981", "pdb1bna"))  # fmt: skip
    add("Structure", "split", "Base pairing: A with T, G with C",
        cls="wide-figure",
        lead="A <b>purine</b> (A or G, two rings) always pairs with a <b>pyrimidine</b> (T or C, "
        "one ring), held by <b>hydrogen bonds</b>.",
        points=[
            "<b>A·T: 2</b> hydrogen bonds. <b>G·C: 3</b> hydrogen bonds, so G·C-rich DNA is "
            "harder to separate.",
            f"Measured here on 1BNA: all {int(m['n_hbonds'])} hydrogen bonds are "
            f"{m['hbond_min']:.2f} to {m['hbond_max']:.2f} Å long.",
            f"Both pairs are equally wide ({m['c1_distance_mean']:.1f} Å from sugar to sugar): "
            "the helix has the same diameter along its whole length.",
        ],
        takeaway="One strand <b>determines</b> the other. This is the basis of heredity.",
        figure={"type": "viewers",
                "items": [
                    {"title": "Adenine · Thymine", "scene": "pair_at",
                     "options": {"colour": "element", "distanceLabels": True, "selectable": False}},
                    {"title": "Guanine · Cytosine", "scene": "pair_gc",
                     "options": {"colour": "element", "distanceLabels": True, "selectable": False}},
                ],
                "caption": "Two real base pairs of 1BNA seen along the helix axis; distances "
                "between the bonded atoms in ångström, measured from the coordinates."},
        sources=_cites("watson1953", "pdb1bna"))  # fmt: skip
    add("Structure", "split", "The double helix: two antiparallel strands",
        lead="The two strands wind around a common axis as a <b>right-handed</b> helix and run "
        "in <b>opposite directions</b>: one 5′→3′, its partner 3′→5′.",
        points=[
            "Outside: the charged sugar-phosphate <b>backbones</b>. Inside: the stacked "
            "<b>base pairs</b>, protected from water.",
            "The two backbones are not opposite each other, which leaves a wide <b>major "
            "groove</b> and a narrow <b>minor groove</b>, where proteins read the sequence.",
        ],
        facts=[
            {"value": f"{m['rise_mean']:.2f} Å", "label": "rise per base pair, measured on 1BNA",
             "source": f"published: {facts['rise']['value']} ({library.cite('watson1953')})",
             "tone": ""},
            {"value": f"{m['bp_per_turn']:.1f}", "label": "base pairs per turn, measured on 1BNA",
             "source": f"in solution: {facts['repeat_solution']['value']} "
             f"({library.cite('wang1979')})", "tone": "good"},
        ],
        figure={"type": "viewer", "scene": "atomic",
                "options": {"colour": "strand", "spin": True, "zoom": 1.0},
                "caption": "All 486 DNA atoms of PDB 1BNA (sequence CGCGAATTCGCG). Click a "
                "nucleotide to see its partner and its hydrogen bonds."},
        sources=_cites("drew1981", "wang1979", "watson1953"))  # fmt: skip
    add("Structure", "text", "The helix in numbers: measured here, compared with the literature",
        html=_table(
            ["Quantity", "Measured on PDB 1BNA by this program", "Published value"],
            [[r["quantity"], r["measured"], r["published"]] for r in biology.helix_rows()],
            numeric=(1,),
        ),
        takeaway="The measurements on real coordinates <b>agree</b> with the published "
        "structure of B-DNA; nothing on these slides is a drawing from memory.",
        sources=_cites("drew1981", "watson1953", "wang1979"))  # fmt: skip
    add("Structure", "split", "Chargaff’s rule, checked on the two genomes of this project",
        lead="If A pairs with T and G with C, a double helix must contain <b>as much A as T and "
        "as much G as C</b>.",
        points=[
            f"<b>phiX174</b>: A {100 * px['A'] / len(phix):.1f}% but T "
            f"{100 * px['T'] / len(phix):.1f}%. The rule fails, and correctly so: inside the "
            "virus this genome is a <b>single strand</b>.",
            f"<b>lambda</b>, a double-stranded genome: on one strand A "
            f"{100 * lm['A'] / len(lam):.1f}% and T {100 * lm['T'] / len(lam):.1f}%, already "
            "close; counted over both strands they are equal by construction.",
        ],
        takeaway="Unequal amounts of A and T are the signature of a <b>single-stranded</b> "
        "genome.",
        figure={"type": "html", "html": diagrams.composition_svg(list(comp.items())),
                "caption": "Counted from the sequences in data/raw (NCBI RefSeq)."},
        sources=_cites("chargaff1950", "refseq_phix174", "refseq_lambda"))  # fmt: skip

    # ---- 2 genomes -------------------------------------------------------------------------------
    add("Genomes", "split", "Genes, genomes and chromosomes",
        cls="wide-figure",
        lead="A <em class='term'>gene</em> is a stretch of DNA that is copied into RNA. The "
        "<em class='term'>genome</em> is all the DNA of an organism; it is carried on one or "
        "more <em class='term'>chromosomes</em>.",
        facts=[_fact("human"), _fact("human_genes", "good"), _fact("chromosomes_human")],
        figure={"type": "html", "html": diagrams.genome_sizes_svg(biology.genome_rows()),
                "caption": "Dark bars: the two genomes used in this project. A human genome is "
                f"about {int(facts['human']['number']) // len(phix):,} times longer than phiX174."},
        sources=_cites("fleischmann1995", "blattner1997", "goffeau1996", "nurk2022", "ihgsc2004"))  # fmt: skip
    add("Genomes", "figure", "Packaging: how a long DNA molecule fits into a nucleus",
        figure={"type": "html", "html": diagrams.packaging_svg(),
                "caption": f"Schematic, not to scale. {facts['nucleosome']['statement']} "
                f"{facts['chromatin_fibre']['statement']} Bacteria and viruses have no "
                "histones: their DNA is compacted by other proteins."},
        takeaway="Eukaryotic DNA is wound around <b>histones</b> into <b>nucleosomes</b> and "
        "folded further; the most compact form is the <b>chromosome</b> seen during cell "
        "division.",
        sources=_cites("kornberg1974", "luger1997", "ou2017", "maeshima2010"))  # fmt: skip
    pair_text = ", ".join(f"{a}/{b}" for a, b, _ in overlap["pairs"][:6])
    add("Genomes", "figure", "A complete genome on one slide: bacteriophage phiX174",
        figure={"type": "html",
                "html": diagrams.gene_map_svg(genes, len(phix), library.phix_annotation()["transcripts"]),
                "caption": f"{overlap['n_genes']} protein-coding genes in {len(phix):,} "
                "nucleotides (NCBI annotation). Each row is one of the three reading frames; "
                f"genes that share DNA lie on different rows ({pair_text} …). ↩ the genome is "
                "a circle: these genes continue at position 1."},
        takeaway=f"The genes add up to <b>{overlap['sum_of_gene_lengths']:,}</b> nucleotides, "
        f"more than the genome has (<b>{len(phix):,}</b>): the same DNA serves <b>more than one "
        "gene</b>, partly in different reading frames. " + str(facts["phix_first"]["value"])
        + ": the first DNA genome ever sequenced.",
        sources=_cites("sanger1977phix", "barrell1976", "refseq_phix174"))  # fmt: skip

    # ---- 3 replication ---------------------------------------------------------------------------
    add("Replication", "figure", "Replication: each strand is the template for a new one",
        figure={"type": "process", "scenario": "replication"},
        sources=_cites("meselson1958", "okazaki1968", "refseq_lambda"))  # fmt: skip
    add("Replication", "text", "Replication in numbers",
        lead="Replication is <b>semiconservative</b>, runs only <b>5′ → 3′</b>, and is fast and "
        "extremely accurate.",
        facts=[_fact("semiconservative", "good"), _fact("fork_speed"), _fact("fidelity", "signal"),
               _fact("okazaki")],
        points=[
            "<b>Why it matters:</b> every cell division must hand a complete, correct copy of "
            "the genome to both daughter cells.",
            "The rare errors that remain are one source of <b>mutations</b> (part 5).",
        ],
        sources=_cites("meselson1958", "okazaki1968", "alberts2022"))  # fmt: skip

    # ---- 4 expression ----------------------------------------------------------------------------
    add("Expression", "split", "The central dogma: DNA → RNA → protein",
        cls="wide-figure",
        lead="A gene is <b>transcribed</b> into RNA; a messenger RNA is <b>translated</b> into "
        "a protein. Proteins do most of the work of a cell.",
        points=[
            "<b>RNA</b> differs from DNA in three ways: the sugar is ribose, the base "
            "<b>uracil (U)</b> replaces thymine, and it is usually a single strand.",
            str(facts["central_dogma"]["statement"]),
        ],
        figure={"type": "html", "html": diagrams.central_dogma_svg(), "caption": ""},
        sources=_cites("crick1970"))  # fmt: skip
    add("Expression", "figure", "Transcription: RNA polymerase copies one strand into RNA",
        figure={"type": "process", "scenario": "transcription"},
        sources=_cites("refseq_phix174", "alberts2022"))  # fmt: skip
    table_css = f"<style>{diagrams.CODON_TABLE_CSS}</style>"
    add("Expression", "split", "The genetic code: 64 codons, 20 amino acids",
        lead="The code is read in <b>triplets</b>. With 4 letters there are 4 × 4 × 4 = "
        f"<b>{len(sequence.CODON_TABLE)}</b> codons.",
        points=[
            f"<b>{len(sequence.CODON_TABLE) - n_stop}</b> codons name an amino acid, "
            f"<b>{n_stop}</b> mean <b>stop</b> ({', '.join(sequence.transcribe(c) for c in sequence.STOP_CODONS)}).",
            "The code is <b>degenerate</b>: most amino acids have several codons, which usually "
            "differ in the third letter.",
            "<b>AUG</b> (methionine) is the start codon. The code is nearly the same in all "
            "living things.",
        ],
        figure={"type": "html", "html": table_css + diagrams.codon_table_html(),
                "caption": "Colours: side-chain class of the amino acid (grey nonpolar, blue "
                "polar, orange acidic, green basic, red stop). Built from NCBI translation "
                "table 1 and compared with Biopython by a test."},
        sources=_cites("crick1961", "nirenberg1961", "ncbi_gencode"))  # fmt: skip
    add("Expression", "figure", "Translation: the ribosome reads codons and builds a protein",
        figure={"type": "process", "scenario": "translation"},
        sources=_cites("nissen2000", "refseq_phix174"))  # fmt: skip

    # ---- 5 mutation and biotechnology ------------------------------------------------------------
    tone = {"silent": "good", "missense": "", "nonsense": "signal", "frameshift": "signal"}
    rows = [
        [r["change"], f"<code>{r['codon']}</code>",
         f'<span class="tag {tone.get(r["effect"], "")}">{r["effect"]}</span>',
         f"{len(r['protein_after'])} of {len(r['protein_before'])} amino acids"
         + ("" if r["protein_after"] != r["protein_before"] else ", unchanged")]
        for r in biology.mutation_rows()
    ]  # fmt: skip
    add("Mutation and biotechnology", "text", "Mutations: one changed letter, four possible outcomes",
        lead="A <em class='term'>mutation</em> is a heritable change of the DNA sequence. Its "
        "effect depends on <b>where</b> it falls and on the <b>reading frame</b>.",
        html=_table(["Change in gene J of phiX174", "Codon", "Effect", "Protein"], rows),
        takeaway="Substitutions change at most one amino acid; an insertion or deletion that is "
        "<b>not a multiple of 3</b> shifts the frame and changes everything after it. "
        "Example in humans: " + str(facts["sickle"]["value"]) + " causes sickle-cell disease.",
        sources=_cites("ingram1957", "refseq_phix174"))  # fmt: skip
    add("Mutation and biotechnology", "text", "Reading, copying and editing DNA",
        lead="Modern biotechnology rests on enzymes that nature uses on DNA.",
        facts=[_fact("sanger_method"), _fact("pcr", "good"), _fact("crispr", "signal")],
        points=[
            f"<b>Sequencing</b>: {facts['sanger_method']['statement']}",
            f"<b>PCR</b>: {facts['pcr']['statement']}",
            f"<b>Genome editing</b>: {facts['crispr']['statement']}",
        ],
        sources=_cites("sanger1977seq", "saiki1988", "jinek2012"))  # fmt: skip

    # ---- 6 sequencing and assembly ---------------------------------------------------------------
    add("Sequencing and assembly", "split", "Sequencing reads DNA only in short pieces",
        lead="In Sanger’s method DNA polymerase copies a template and stops at a random "
        "position of one chosen letter. Sorting the copies by length <b>spells the sequence</b>.",
        points=[
            "One reaction reads only a short stretch, never a whole genome.",
            "Modern machines read millions of pieces in parallel, but each <b>read</b> is "
            "still short, and a few letters are wrong.",
        ],
        takeaway="A genome must be <b>reconstructed</b> from many short, overlapping, slightly "
        "wrong reads.",
        figure={"type": "html", "html": diagrams.sanger_gel_svg(gel_letters, gel_first),
                "caption": "Idealised gel computed from the real sequence at the start of "
                "gene J: one band for every stopped copy."},
        sources=_cites("sanger1977seq"))  # fmt: skip
    add("Sequencing and assembly", "figure", "Shotgun sequencing: random reads that overlap",
        figure={"type": "html",
                "html": diagrams.shotgun_svg(shotgun["starts"], shotgun["lengths"], shotgun["length"]),
                "caption": f"A real simulation: {shotgun['n']} reads of {shotgun['read_length']} "
                f"letters from the first {shotgun['length']} letters of phiX174, coverage "
                f"{shotgun['coverage']:.0f}×. The positions are known here; a sequencer does "
                "not report them."},
        takeaway="<b>Coverage</b> = reads × read length ÷ genome length. Enough coverage makes "
        "neighbouring reads <b>overlap</b>, and overlaps are the only clue to their order.",
        sources=_cites("fleischmann1995", "lander1988"))  # fmt: skip
    add("Sequencing and assembly", "split", "Assembly is a travelling salesman problem",
        lead="Put the reads in an order and glue each to the next. The letters saved by the "
        "overlaps are the <b>fitness F</b>; the best order saves the most.",
        html='<div class="formula">F(π) = Σ ov(π<sub>k</sub>, π<sub>k+1</sub>) &nbsp;→&nbsp; max'
        "<small>ov(i, j): how many letters at the end of read i equal the start of read j"
        "</small></div>",
        points=[
            f"Toy example: reads {', '.join(toy['reads'])}. Order <b>{toy['order']}</b> saves "
            f"<b>{toy['f']}</b> letters: {toy['total']} letters shrink to {toy['l']} "
            f"(<code>{toy['merged']}</code>).",
            "Visiting every read once along the best path is the <b>asymmetric travelling "
            "salesman problem</b>: overlap from A to B differs from B to A.",
            "The number of orders is n!, so exhaustive search is impossible beyond a few reads.",
        ],
        figure={"type": "image", "src": toy_uri, "alt": "Overlap graph of the toy example",
                "caption": "Arrows: overlaps; red: the best order."})  # fmt: skip
    add("Sequencing and assembly", "text", "Four ways to search for the best order",
        html=_table(
            ["Method", "Idea", "Role in this project"],
            [
                ["Greedy", "Always glue the two pieces with the largest overlap", "fast baseline"],
                ["Tabu search", "Improve step by step and forbid recent moves, to leave local "
                 "optima", "classic local search"],
                ["LKH-3", "The strongest known heuristic for travelling salesman problems",
                 "reference solver"],
                ["ALNS + bandit", "Destroy part of the order and repair it with one of six "
                 "operator pairs; a <b>multi-armed bandit</b> learns which pair to use",
                 "the method under test"],
            ],
        ),
        takeaway="Research question, fixed before the experiments: does a <b>bandit</b> choose "
        "operators better than <b>random choice</b> under the same time budget?",
        sources=_cites("glover1989", "helsgaun2000", "ropke2006", "auer2002"))  # fmt: skip
    add("Sequencing and assembly", "split", "Result: an honest negative answer",
        cls="wide-figure",
        lead=f"<b>{res['n_sets']}</b> read sets ({res['n_settings']} settings × {res['seeds']} "
        f"seeds), every solver with the same <b>{res['budget']:.0f} seconds</b>.",
        points=[
            f"<b>The pre-registered claim is not supported</b>: UCB1 is not different from "
            f"random choice ({res['ucb_p']}, A12 = {res['ucb_a12']:.2f}); epsilon-greedy is "
            f"slightly worse (A12 = {res['eps_a12']:.2f}).",
            "LKH-3 was <b>never beaten</b>; ALNS beats greedy and tabu search."
            if res["lkh_never_beaten"] else "ALNS beats greedy and tabu search.",
            f"<b>Post hoc, exploratory only</b>: UCB1 with random tie-breaking is slightly "
            f"better than random (A12 = {res['post_a12']:.2f}, {res['post_p']}). It was added "
            "after seeing the results and is not a confirmation.",
        ],
        figure={"type": "html",
                "html": _table(["Solver", "Mean gap to best", "Contigs", "Identity"], res["rows"],
                               numeric=(1, 2, 3)),
                "caption": "Gap: how far the fitness is below the best solution found on the "
                "same read set (lower is better). Contigs and identity are measured against "
                "the real genome. Source: results/phase5."})  # fmt: skip

    # ---- end -------------------------------------------------------------------------------------
    add("Summary", "summary", "What to remember", kicker="Summary", points=[
        "DNA is a <b>double helix of two antiparallel strands</b>; A pairs with T (2 hydrogen "
        "bonds), G with C (3).",
        f"The helix measured here on real atoms: <b>{m['rise_mean']:.2f} Å</b> per base pair, "
        f"<b>{m['bp_per_turn']:.1f}</b> base pairs per turn, right-handed.",
        "<b>Replication</b> is semiconservative and runs 5′→3′: one strand continuously, the "
        "other in Okazaki fragments.",
        "<b>Transcription</b> copies a gene into RNA; <b>translation</b> reads codons of three "
        f"letters: gene J of phiX174 gives a protein of {len(protein_j)} amino acids.",
        "A <b>mutation</b> can be silent, change one amino acid, stop the protein early or "
        "shift the whole reading frame.",
        "Sequencing gives short reads; <b>assembly</b> orders them by their overlaps, which is "
        "an asymmetric <b>travelling salesman problem</b>.",
    ])  # fmt: skip
    for part, page in enumerate(_references_html(), 1):
        add("References", "references", f"References ({part} of 2)", kicker="Sources", html=page)

    first_of: dict[str, int] = {}
    for i, slide in enumerate(slides):
        first_of.setdefault(slide["section"], i)
    slides[1]["items"] = [
        {"name": name, "what": what, "slide": first_of[name]} for name, what in SECTIONS
    ]
    return {
        "title": TITLE,
        "slides": slides,
        "scenes": biology.scenes(),
        "scenarios": biology.scenarios(),
    }


def html(author: str = "", *, embedded: bool = False) -> str:
    return studio.lecture_html(build(author), embedded=embedded)


def write(path: Path, author: str = "") -> Path:
    """Write the lecture as one HTML file that opens by double-click, without the application."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html(author), encoding="utf-8", newline="\n")
    return path


__all__ = ["FIGURES", "SECTIONS", "TITLE", "build", "html", "write"]
