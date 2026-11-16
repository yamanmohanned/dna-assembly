"""Molecular biology layer: sequences, the genetic code, mutations, the 1BNA structure, the scenes.

Expected values are justified by theory or by an independent implementation (Biopython), never
copied from a first run: the codon table is compared with Biopython's, the helix measurements with
the published ranges for B-DNA, the genes with the definition of a coding sequence.
"""

from __future__ import annotations

import json
import math
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq
from hypothesis import given
from hypothesis import strategies as st

from genomeasm.appkit import lab
from genomeasm.molbio import diagrams, library, model, processes, quiz, sequence
from genomeasm.molbio import structure as structure_mod

PHIX = "NC_001422.1"
LAMBDA = "NC_001416.1"
dna = st.text(alphabet="ACGT", min_size=1, max_size=200)


# ---- sequences and the genetic code --------------------------------------------------------------


def test_codon_table_equals_the_standard_table_of_biopython() -> None:
    standard = CodonTable.unambiguous_dna_by_id[1]
    assert len(sequence.CODON_TABLE) == 64
    for codon, amino in sequence.CODON_TABLE.items():
        expected = "*" if codon in standard.stop_codons else standard.forward_table[codon]
        assert amino == expected, codon
    assert set(sequence.STOP_CODONS) == set(standard.stop_codons) == {"TAA", "TAG", "TGA"}
    # 20 amino acids and the stop signal, every one with its three-letter code
    assert set(sequence.CODON_TABLE.values()) == set(sequence.AMINO_ACIDS)
    assert len(sequence.AMINO_ACIDS) == 21


def test_bacterial_table_has_the_same_meanings_and_more_start_codons() -> None:
    bacterial = CodonTable.unambiguous_dna_by_id[11]
    assert bacterial.forward_table == CodonTable.unambiguous_dna_by_id[1].forward_table
    assert set(sequence.BACTERIAL_STARTS) == set(bacterial.start_codons)


@given(dna)
def test_reverse_complement_is_its_own_inverse_and_matches_biopython(s: str) -> None:
    assert sequence.reverse_complement(sequence.reverse_complement(s)) == s
    assert sequence.reverse_complement(s) == str(Seq(s).reverse_complement())
    assert sequence.complement(s) == str(Seq(s).complement())
    assert sequence.transcribe(s) == str(Seq(s).transcribe())


@given(dna)
def test_translation_matches_biopython(s: str) -> None:
    whole = s[: len(s) - len(s) % 3]
    assert sequence.translate(whole, to_stop=False) == str(Seq(whole).translate())
    assert sequence.translate(whole) == str(Seq(whole).translate(to_stop=True))


@given(dna)
def test_hydrogen_bonds_are_two_per_at_and_three_per_gc(s: str) -> None:
    at, gc = s.count("A") + s.count("T"), s.count("G") + s.count("C")
    assert sequence.hydrogen_bonds(s) == 2 * at + 3 * gc
    assert sequence.gc_fraction(s) == pytest.approx(gc / len(s))


def test_anticodon_pairs_antiparallel_with_its_codon() -> None:
    assert sequence.anticodon("AUG") == "CAU"  # 5'-CAU-3' pairs with 5'-AUG-3'
    assert sequence.clean_dna(" augGc x") == "ATGGC"
    with pytest.raises(ValueError):
        sequence.complement("ACGU")


def test_open_reading_frames_are_found_in_every_frame() -> None:
    #        frame 0: ATG AAA CCC TAA          frame 1 (from letter 14): ATG GGG TGA
    text = "ATGAAACCCTAA" + "C" + "ATGGGGTGA"
    orfs = sequence.find_orfs(text, min_codons=2)
    assert [(o.start, o.end, o.protein, o.frame) for o in orfs] == [
        (0, 12, "MKP", 0),
        (13, 22, "MG", 1),
    ]
    assert sequence.find_orfs(text, min_codons=3) == [orfs[0]]


# ---- mutations -----------------------------------------------------------------------------------

CDS = "ATGTCTAAAGGTTAA"  # Met Ser Lys Gly stop


@pytest.mark.parametrize(
    ("kind", "position", "base", "effect", "protein"),
    [
        ("substitution", 5, "C", "silent", "MSKG"),  # TCT -> TCC, both serine
        ("substitution", 6, "G", "missense", "MSEG"),  # AAA -> GAA, lysine -> glutamate
        ("substitution", 6, "T", "nonsense", "MS"),  # AAA -> TAA, a stop codon
        ("substitution", 1, "A", "start lost", "KSKG"),  # ATG -> AAG, not a start codon
        ("substitution", 12, "C", "stop lost", "MSKGQ"),  # TAA -> CAA
        ("deletion", 4, "", "frameshift", "MLKV"),  # ATG TTA AAG GTT AA: no stop codon left
        ("insertion", 6, "GGC", "in-frame insertion", "MSGKG"),
    ],
)
def test_mutations_are_classified_by_the_textbook_definitions(
    kind: str, position: int, base: str, effect: str, protein: str
) -> None:
    mutated = sequence.apply_mutation(CDS, kind, position, base)
    result = sequence.classify_mutation(CDS, mutated)
    assert result.effect == effect
    assert result.protein_before == "MSKG"
    assert result.protein_after == protein == sequence.translate(mutated)


def test_a_mutation_must_change_something_and_stay_inside_the_sequence() -> None:
    with pytest.raises(ValueError, match="nothing is mutated"):
        sequence.apply_mutation(CDS, "substitution", 0, "A")
    with pytest.raises(ValueError, match="outside"):
        sequence.apply_mutation(CDS, "deletion", len(CDS))
    with pytest.raises(ValueError):
        sequence.apply_mutation(CDS, "inversion", 0, "A")
    assert sequence.classify_mutation(CDS, CDS).effect == "no change"
    assert sequence.substitution_class("A", "G") == "transition"
    assert sequence.substitution_class("A", "T") == "transversion"


@given(st.integers(0, len(CDS) - 1), st.sampled_from("ACGT"))
def test_a_deletion_or_single_insertion_always_shifts_the_frame(position: int, base: str) -> None:
    deleted = sequence.apply_mutation(CDS, "deletion", position)
    inserted = sequence.apply_mutation(CDS, "insertion", position, base)
    assert sequence.classify_mutation(CDS, deleted).effect == "frameshift"
    assert sequence.classify_mutation(CDS, inserted).effect == "frameshift"


# ---- the annotated genes of phiX174 --------------------------------------------------------------


def test_every_annotated_gene_is_a_clean_coding_sequence_of_the_real_genome() -> None:
    genome = lab.genome_sequence(PHIX)
    genes = library.phix_genes()
    assert len(genes) == 11 and len({g.name for g in genes}) == 11
    for gene in genes:
        cds = gene.sequence(genome)
        assert all(1 <= a <= b <= len(genome) for a, b in gene.parts), gene.name
        assert len(cds) == gene.length and len(cds) % 3 == 0, gene.name
        assert cds[:3] in sequence.BACTERIAL_STARTS, gene.name
        protein = sequence.translate(cds, to_stop=False)
        assert protein.count("*") == 1 and protein.endswith("*"), gene.name
        assert len(sequence.translate(cds)) == len(cds) // 3 - 1
        assert cds == str(Seq(cds)) and protein == str(Seq(cds).translate(table=11))


def test_genes_of_phix174_overlap_as_reported_by_sanger() -> None:
    genome = lab.genome_sequence(PHIX)
    genes = library.phix_genes()
    pairs = {(a, b) for a, b, _ in library.overlapping_pairs(genes)}
    # "two pairs of genes are coded by the same region of DNA" (Sanger et al. 1977): B in A, E in D
    assert ("A", "B") in pairs and ("D", "E") in pairs
    assert library.gene("B").frame() != library.gene("A").frame()
    assert library.gene("E").frame() != library.gene("D").frame()
    # so the genes add up to more letters than the genome has
    assert sum(g.length for g in genes) > len(genome)


# ---- the experimental structure ------------------------------------------------------------------


@pytest.fixture(scope="module")
def bna() -> structure_mod.Structure:
    return library.structure_1bna()


def test_pdb_file_is_read_completely(bna: structure_mod.Structure) -> None:
    text = library.PDB_1BNA.read_text(encoding="ascii")
    assert len(bna.atoms) == sum(1 for line in text.splitlines() if line.startswith("ATOM"))
    assert bna.pdb_id == "1BNA" and bna.resolution == pytest.approx(1.9)
    assert "X-RAY" in bna.experiment and bna.doi == "10.1073/pnas.78.4.2179"
    assert bna.chains() == ["A", "B"] and len(bna.residues) == 24
    # the molecule named in the file header: d(CGCGAATTCGCG), a self-complementary sequence
    assert "CP*GP*CP*GP*AP*AP*TP*TP*CP*GP*CP*G" in text
    for chain in bna.chains():
        assert bna.sequence(chain) == "CGCGAATTCGCG"
    assert sequence.reverse_complement("CGCGAATTCGCG") == "CGCGAATTCGCG"
    assert {a.element for a in bna.atoms} == {"C", "N", "O", "P"}  # no hydrogens at 1.9 angstrom


def test_covalent_bonds_have_chemical_lengths(bna: structure_mod.Structure) -> None:
    bonds = structure_mod.covalent_bonds(bna)
    lengths = [float(np.linalg.norm(bna.coords[i] - bna.coords[j])) for i, j in bonds]
    assert min(lengths) > 1.1 and max(lengths) < 1.95  # C=O about 1.2, P-O about 1.6 angstrom
    bonded = {i for bond in bonds for i in bond}
    assert bonded == set(range(len(bna.atoms)))  # no atom floats free
    links = [
        (i, j) for i, j in bonds
        if bna.atoms[i].resseq != bna.atoms[j].resseq
    ]  # fmt: skip
    assert len(links) == 22  # 12 nucleotides per chain are joined by 11 phosphodiester links
    assert all({bna.atoms[i].name, bna.atoms[j].name} == {"O3'", "P"} for i, j in links)


def test_watson_crick_pairs_and_hydrogen_bonds_are_found_from_the_coordinates(
    bna: structure_mod.Structure,
) -> None:
    pairs = structure_mod.base_pairs(bna)
    assert len(pairs) == 12
    for pair in pairs:
        purine, pyrimidine = bna.residues[pair.purine], bna.residues[pair.pyrimidine]
        assert sequence.COMPLEMENT[purine.base] == pyrimidine.base
        assert purine.chain != pyrimidine.chain
        assert purine.resseq + pyrimidine.resseq == 25  # antiparallel: residue i pairs with 25 - i
        assert len(pair.hbonds) == sequence.HYDROGEN_BONDS[purine.base]  # 2 for A-T, 3 for G-C
        for bond in pair.hbonds:
            assert 2.5 < bond.distance < 3.5  # N-H...O and N-H...N bonds are about 2.7 to 3.1
    assert sum(len(p.hbonds) for p in pairs) == 2 * 4 + 3 * 8  # four A-T and eight G-C pairs


def test_helix_measurements_agree_with_published_b_dna(bna: structure_mod.Structure) -> None:
    geometry = structure_mod.helix_geometry(bna)
    summary = geometry.summary()
    assert geometry.n_pairs == 12 and len(geometry.rises) == len(geometry.twists) == 11
    assert all(t > 0 for t in geometry.twists)  # positive twist: a right-handed helix
    assert summary["rise_mean"] == pytest.approx(3.4, abs=0.15)  # Watson and Crick: 3.4 angstrom
    # 10.0 base pairs per turn in fibres, 10.4 in solution, 9.8 in the centre of this crystal
    assert 9.6 < summary["bp_per_turn"] < 10.6
    assert summary["pitch"] == pytest.approx(34.0, abs=2.0)
    assert 8.0 < summary["p_radius_mean"] < 10.0  # phosphorus about 1 nm from the axis
    assert summary["c1_distance_mean"] == pytest.approx(10.5, abs=0.4)
    assert math.isclose(float(np.linalg.norm(geometry.axis)), 1.0, abs_tol=1e-9)
    c1 = geometry.c1_distances  # A-T and G-C pairs are equally wide
    assert max(c1) - min(c1) < 1.0


def test_backbone_template_puts_phosphate_outside_and_the_base_inside(
    bna: structure_mod.Structure,
) -> None:
    template = structure_mod.backbone_template(bna)
    radius = {key: math.hypot(value[0], value[1]) for key, value in template.items()}
    assert radius["P"] > radius["sugar"] > radius["C1'"] > radius["edge"]
    assert radius["edge"] < 2.0 and radius["P"] > 8.0


# ---- scenes for the 3D viewer --------------------------------------------------------------------


def _phi(scene: dict, strand: int, element: str) -> list[float]:
    atoms = scene["atoms"]
    xyz = np.array(atoms["xyz"]).reshape(-1, 3)
    chosen = [
        xyz[i]
        for i in range(len(xyz))
        if atoms["el"][i] == element and scene["residues"][atoms["res"][i]]["strand"] == strand
    ]
    chosen.sort(key=lambda p: p[1])
    return [math.degrees(math.atan2(p[0], p[2])) for p in chosen]


def _turns_right(angles: list[float]) -> bool:
    steps = [(b - a + 180) % 360 - 180 for a, b in zip(angles, angles[1:], strict=False)]
    return all(step > 0 for step in steps)


def test_atomic_scene_keeps_every_atom_and_is_right_handed(bna: structure_mod.Structure) -> None:
    scene = model.atomic_scene(bna)
    n = len(bna.atoms)
    assert len(scene["atoms"]["el"]) == n and len(scene["atoms"]["xyz"]) == 3 * n
    assert len(scene["hbonds"]) == 32 and len(scene["bonds"]) % 2 == 0
    assert {label["text"] for label in scene["labels"]} >= {
        "5′",
        "3′",
        "major groove",
        "minor groove",
    }
    # going up the axis (+y), the phosphates of one strand turn by the right-hand rule
    assert _turns_right(_phi(scene, 0, "P"))
    json.dumps(scene)  # plain data, ready for the browser


def test_subset_scene_shows_one_pair_flat_with_its_measured_bonds(
    bna: structure_mod.Structure,
) -> None:
    scene = model.atomic_scene(bna, [4, 19], orient="pair")
    xyz = np.array(scene["atoms"]["xyz"]).reshape(-1, 3)
    base_atoms = [i for i, c in enumerate(scene["atoms"]["comp"]) if c == 2]
    assert np.ptp(xyz[base_atoms, 2]) < 1.5  # the bases lie in the plane of the screen
    assert len(scene["hbonds"]) == 2 and not scene["labels"]
    assert {h["text"] for h in scene["hbonds"]} == {"A N6 ··· O4 T", "A N1 ··· N3 T"}
    assert scene["extent"][0] > scene["extent"][1] > scene["extent"][2]


def test_schematic_scene_follows_the_sequence_and_the_canonical_geometry() -> None:
    piece = lab.genome_sequence(PHIX)[:21]
    scene = model.schematic_scene(piece, library.backbone_template())
    residues = scene["residues"]
    assert [r["base"] for r in residues[:21]] == list(piece)
    assert [r["base"] for r in residues[21:]] == list(sequence.complement(piece))
    assert len(scene["hbonds"]) == sequence.hydrogen_bonds(piece)
    assert _turns_right(_phi(scene, 0, "P")) and _turns_right(_phi(scene, 1, "P"))
    xyz = np.array(scene["atoms"]["xyz"]).reshape(-1, 3)
    heights = sorted(
        xyz[i][1] for i in range(len(xyz))
        if scene["atoms"]["el"][i] == "P" and residues[scene["atoms"]["res"][i]]["strand"] == 0
    )  # fmt: skip
    assert np.allclose(np.diff(heights), model.SCHEMATIC_RISE, atol=1e-3)
    steps = np.diff(np.unwrap(np.radians(_phi(scene, 0, "P"))))
    assert np.allclose(np.degrees(steps), 360 / model.SCHEMATIC_BP_PER_TURN, atol=0.01)
    with pytest.raises(ValueError):
        model.schematic_scene("A", library.backbone_template())


# ---- references, facts, quiz ---------------------------------------------------------------------


def test_references_are_complete_and_facts_point_to_them() -> None:
    refs = library.references()
    assert len(refs) >= 40
    dois = [r["doi"] for r in refs.values() if r["doi"]]
    assert len(dois) == len(set(dois))
    for key, ref in refs.items():
        assert ref["authors"] and ref["title"] and ref["journal"], key
        assert ref["check"] in library.CHECK_LABELS and ref["topic"] in library.TOPIC_ORDER, key
        if ref["check"] == "crossref":
            assert ref["doi"] and ref["volume"] and 1940 <= int(ref["year"]) <= 2026, key
        family = str(ref["authors"][0]).split(" ")[0]
        assert family in library.reference_line(ref), key
        assert ref.get("short") or (
            family in library.cite(key) and str(ref["year"]) in library.cite(key)
        ), key
    for key, fact in library.facts().items():
        assert fact["source"] in refs, key
        assert fact["check"] in ("abstract", "corroborated", "textbook"), key
        assert fact["value"] and fact["label"] and fact["statement"], key
    assert library.cite("watson1953") == "Watson & Crick 1953"
    assert library.cite("drew1981") == "Drew et al. 1981"
    assert library.cite("pdb1bna") == "PDB 1BNA"


def test_quoted_genome_sizes_carry_the_number_of_their_source() -> None:
    facts = library.facts()
    assert facts["ecoli"]["number"] == 4_639_221 and "4,639,221" in facts["ecoli"]["value"]
    assert facts["hinfluenzae"]["number"] == 1_830_137
    assert facts["human"]["number"] == 3_055_000_000 and "3.055" in facts["human"]["value"]
    assert len(lab.genome_sequence(PHIX)) == 5386 and len(lab.genome_sequence(LAMBDA)) == 48502


def test_every_quiz_question_has_one_valid_computed_answer() -> None:
    questions = quiz.questions()
    assert len(questions) >= 10
    assert {q.topic for q in questions} == {"structure", "replication", "expression", "mutation"}
    for q in questions:
        assert 0 <= q.answer < len(q.options) and len(set(q.options)) == len(q.options) >= 3
        assert q.why and q.text
    by_text = {q.text: q for q in questions}
    partner = next(q for t, q in by_text.items() if "other strand" in t)
    expected = str(Seq(quiz.EXAMPLE).reverse_complement())
    assert partner.options[partner.answer] == f"5′-{expected}-3′"
    protein = next(q for t, q in by_text.items() if t.startswith("Translate"))
    assert protein.options[protein.answer] == "Met-Ser-Lys-Gly"
    assert lab.genome_sequence(PHIX)[847:859] == quiz.EXAMPLE  # the real start of gene J


# ---- animation scenarios -------------------------------------------------------------------------


@pytest.fixture(scope="module")
def scenarios() -> dict[str, dict]:
    return processes.scenarios(lab.genome_sequence(LAMBDA), lab.genome_sequence(PHIX))


def test_every_scenario_is_a_consistent_timeline(scenarios: dict[str, dict]) -> None:
    assert set(scenarios) == {"replication", "transcription", "translation"}
    for name, scenario in scenarios.items():
        steps = scenario["steps"]
        keys = set(steps[0]["to"])
        assert len(steps) >= 6 and scenario["note"].startswith("Schematic"), name
        for step in steps:
            assert set(step["to"]) == keys and step["dur"] > 0 and step["title"], name
            for key, (start, end) in step["win"].items():
                assert key in keys and 0 <= start < end <= 1, name
        json.dumps(scenario)


def test_replication_scenario_copies_both_strands_of_a_real_stretch(
    scenarios: dict[str, dict],
) -> None:
    data = scenarios["replication"]["data"]
    start = data["start"]
    top = lab.genome_sequence(LAMBDA)[start - 1 : start - 1 + processes.REPLICATION_LENGTH]
    assert data["top"] == top and data["bottom"] == sequence.complement(top)
    covered = [i for a, b in data["fragments"] for i in range(a, b + 1)]
    assert covered == list(range(len(top)))  # the Okazaki fragments tile the lagging strand
    steps = scenarios["replication"]["steps"]
    for step in steps:  # a strand is never copied ahead of the fork
        assert step["to"]["lead"] <= step["to"]["fork"]
    last = steps[-1]["to"]
    assert last["lead"] == len(top) and last["ligate"] == 1 and last["apart"] == 1
    assert all(last[f"f{k}"] == b - a + 1 for k, (a, b) in enumerate(data["fragments"]))


def test_transcription_and_translation_use_the_real_gene_j(scenarios: dict[str, dict]) -> None:
    genome = lab.genome_sequence(PHIX)
    gene = library.gene("J")
    cds = gene.sequence(genome)
    t = scenarios["transcription"]["data"]
    flank = t["flank"]
    assert t["coding"][flank:].startswith(cds[: processes.TRANSCRIPT_LENGTH])
    assert t["rna"] == str(Seq(cds[: processes.TRANSCRIPT_LENGTH]).transcribe())
    assert t["template"] == sequence.complement(t["coding"])
    rows = scenarios["translation"]["data"]["codons"]
    protein = str(Seq(cds).translate(table=11, to_stop=True))
    assert scenarios["translation"]["data"]["protein"] == protein and len(protein) == 38
    assert "".join(r["one"] for r in rows) == protein + "*"
    assert "".join(r["codon"] for r in rows) == str(Seq(cds).transcribe())
    for row in rows:  # the anticodon is written under the codon: position by position its partner
        pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}
        assert row["anticodon"] == "".join(pairs[c] for c in row["codon"])
    last = scenarios["translation"]["steps"][-1]["to"]
    assert last["cycle"] == len(protein) - 1 and last["free"] == 1


# ---- diagrams ------------------------------------------------------------------------------------


def test_diagrams_are_well_formed_and_drawn_from_their_data() -> None:
    genome = lab.genome_sequence(PHIX)
    genes = library.phix_genes()
    drawings = {
        "nucleotide": diagrams.nucleotide_svg(),
        "dogma": diagrams.central_dogma_svg(),
        "packaging": diagrams.packaging_svg(),
        "genes": diagrams.gene_map_svg(
            genes, len(genome), library.phix_annotation()["transcripts"]
        ),
        "sizes": diagrams.genome_sizes_svg([("phiX174", len(genome), "this project")]),
        "bases": diagrams.composition_svg([("phiX174", sequence.base_composition(genome))]),
        "gel": diagrams.sanger_gel_svg(genome[847:865], 848),
        "shotgun": diagrams.shotgun_svg([0, 30, 70], [60, 60, 30], 120),
    }
    for name, svg in drawings.items():
        root = ET.fromstring(svg)
        assert root.tag.endswith("svg") and root.get("role") == "img" and root.get("aria-label"), (
            name
        )
    for gene in genes:
        assert f">{gene.name}<" in drawings["genes"], gene.name
    assert drawings["gel"].count('height="8"') == 18  # one band per stopped copy
    assert "20 letters are covered by no read" in drawings["shotgun"]  # letters 100 to 119
    assert f"{len(genome):,} bp" in drawings["sizes"]
    table = diagrams.codon_table_html()
    assert table.count('<td class="codon') == 64  # one cell per codon
    assert table.count("<b>Stop</b>") == 3 and "<code>AUG</code> <b>Met</b>" in table
