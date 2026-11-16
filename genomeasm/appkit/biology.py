"""What the molecular-biology pages and the lecture show, computed once from the project data.

No Streamlit import and no hand-typed result: the 3D scenes come from PDB 1BNA and the phiX174
sequence, the tables from measurements and from ``docs/molbio_facts.yaml``, the mutation examples
from ``sequence.classify_mutation`` applied to a real gene.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

from genomeasm.appkit import lab
from genomeasm.molbio import library, model, processes, sequence

PHIX = "NC_001422.1"
LAMBDA = "NC_001416.1"
SCHEMATIC_LETTERS = 32  # about three turns of the helix
SANGER_START = 848  # the start of gene J
SANGER_LETTERS = 18
# residues of 1BNA shown alone (0-based index in the structure): chain A residue 5 is an adenine
NUCLEOTIDE_RESIDUE = 4
PAIR_AT = (4, 19)  # dA5 of chain A with dT20 of chain B
PAIR_GC = (3, 20)  # dG4 of chain A with dC21 of chain B
# point mutations of gene J used as examples: (kind, 0-based position in the gene, new base)
MUTATION_EXAMPLES = (
    ("substitution", 5, "C"),
    ("substitution", 6, "G"),
    ("substitution", 6, "T"),
    ("deletion", 4, ""),
    ("insertion", 6, "GGC"),
)


def genome(accession: str) -> str:
    return lab.genome_sequence(accession)


@lru_cache(maxsize=1)
def scenes() -> dict[str, dict[str, Any]]:
    """The 3D scenes: real atoms of 1BNA and the schematic helix of the phiX174 sequence."""
    structure = library.structure_1bna()
    return {
        "atomic": model.atomic_scene(structure),
        "nucleotide": model.atomic_scene(structure, [NUCLEOTIDE_RESIDUE], orient="pair"),
        "pair_at": model.atomic_scene(structure, list(PAIR_AT), orient="pair"),
        "pair_gc": model.atomic_scene(structure, list(PAIR_GC), orient="pair"),
        "schematic": schematic(PHIX, 1),
    }


def schematic(accession: str, first: int, letters: int = SCHEMATIC_LETTERS) -> dict[str, Any]:
    """Schematic helix of a stretch of a genome (``first`` is a 1-based position)."""
    piece = genome(accession)[first - 1 : first - 1 + letters]
    return model.schematic_scene(piece, library.backbone_template(), first_position=first)


@lru_cache(maxsize=1)
def scenarios() -> dict[str, dict[str, Any]]:
    return processes.scenarios(genome(LAMBDA), genome(PHIX))


def helix_rows() -> list[dict[str, str]]:
    """The double helix in numbers: measured here on 1BNA, next to the published values."""
    m = library.measured_helix()
    facts = library.facts()

    def published(key: str) -> str:
        return f"{facts[key]['value']} ({library.cite(str(facts[key]['source']))})"

    return [
        {
            "quantity": "Rise per base pair",
            "measured": f"{m['rise_mean']:.2f} ± {m['rise_sd']:.2f} Å",
            "published": published("rise"),
        },
        {
            "quantity": "Rotation per base pair (twist)",
            "measured": f"{m['twist_mean']:.1f} ± {m['twist_sd']:.1f}°, right-handed",
            "published": published("repeat_crystal"),
        },
        {
            "quantity": "Base pairs per turn",
            "measured": f"{m['bp_per_turn']:.1f}",
            "published": f"{facts['repeat_fibre']['value']} in fibres, "
            f"{facts['repeat_solution']['value']} in solution ({library.cite('wang1979')})",
        },
        {
            "quantity": "Phosphorus atoms from the axis",
            "measured": f"{m['p_radius_mean']:.1f} Å (so about {m['p_diameter'] / 10:.1f} nm across)",
            "published": published("diameter"),
        },
        {
            "quantity": "Watson-Crick hydrogen bonds",
            "measured": f"{int(m['n_hbonds'])} in {int(m['n_pairs'])} pairs, "
            f"{m['hbond_min']:.2f} to {m['hbond_max']:.2f} Å (mean {m['hbond_mean']:.2f} Å)",
            "published": "2 per A·T pair, 3 per G·C pair",
        },
        {
            "quantity": "Sugar to sugar across a pair (C1′ to C1′)",
            "measured": f"{m['c1_distance_mean']:.1f} Å for A·T and G·C alike",
            "published": "one purine + one pyrimidine: every pair has the same width",
        },
    ]


def genome_rows() -> list[tuple[str, int, str]]:
    """Genome sizes for the comparison chart: (name, base pairs, where the number comes from)."""
    facts = library.facts()
    rows = [
        ("phiX174 (virus)", len(genome(PHIX)), "this project, NCBI NC_001422.1"),
        ("lambda (virus)", len(genome(LAMBDA)), "this project, NCBI NC_001416.1"),
    ]
    for key, name in (
        ("hinfluenzae", "Haemophilus influenzae"),
        ("ecoli", "Escherichia coli K-12"),
        ("yeast", "Yeast (S. cerevisiae)"),
        ("human", "Human (one set)"),
    ):
        rows.append((name, int(facts[key]["number"]), library.cite(str(facts[key]["source"]))))
    return rows


def composition_rows() -> list[tuple[str, dict[str, int]]]:
    return [
        ("phiX174: single-stranded in the virus", sequence.base_composition(genome(PHIX))),
        ("lambda: double-stranded, one strand counted", sequence.base_composition(genome(LAMBDA))),
    ]


def gene_sequence(name: str = "J") -> str:
    return library.gene(name).sequence(genome(PHIX))


def mutation_rows(gene_name: str = "J") -> list[dict[str, Any]]:
    """The example mutations of a real gene with their computed effect on the protein."""
    cds = gene_sequence(gene_name)
    rows = []
    for kind, position, base in MUTATION_EXAMPLES:
        mutated = sequence.apply_mutation(cds, kind, position, base)
        effect = sequence.classify_mutation(cds, mutated)
        if kind == "substitution":
            change = f"{cds[position]} → {base} at letter {position + 1}"
        elif kind == "deletion":
            change = f"{cds[position]} deleted at letter {position + 1}"
        else:
            change = f"{base} inserted before letter {position + 1}"
        rows.append(
            {
                "kind": kind,
                "change": change,
                "effect": effect.effect,
                "codon": f"{effect.codon_before} → {effect.codon_after}",
                "protein_before": effect.protein_before,
                "protein_after": effect.protein_after,
                "explanation": effect.explanation,
            }
        )
    return rows


def sanger_window() -> tuple[str, int]:
    """A real stretch of phiX174 for the sequencing-gel drawing: (letters, 1-based position)."""
    return genome(PHIX)[SANGER_START - 1 : SANGER_START - 1 + SANGER_LETTERS], SANGER_START


@lru_cache(maxsize=1)
def shotgun_example() -> dict[str, Any]:
    """A small real shotgun experiment: reads simulated from a window of phiX174."""
    params = lab.InstanceParams(
        genome=PHIX, start=0, length=600, coverage=6.0, read_length=60, error=0.0, seed=1
    )
    reads = lab.simulate(params)
    return {
        "starts": [int(s) for s in reads.starts],
        "lengths": [int(s) for s in reads.spans],
        "length": params.length,
        "n": reads.n_reads,
        "coverage": params.coverage,
        "read_length": params.read_length,
    }


def overlap_summary() -> dict[str, Any]:
    """Overlapping genes of phiX174, counted from the annotation."""
    genes = library.phix_genes()
    pairs = library.overlapping_pairs(genes)
    coding = {p for g in genes for a, b in g.parts for p in range(a, b + 1)}
    return {
        "n_genes": len(genes),
        "pairs": pairs,
        "coding_letters": len(coding),
        "genome_length": len(genome(PHIX)),
        "sum_of_gene_lengths": sum(g.length for g in genes),
    }
