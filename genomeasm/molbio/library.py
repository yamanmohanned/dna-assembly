"""References, quoted facts and gene annotation of the molecular-biology pages (read from YAML).

Nothing here is typed twice: the pages, the lecture and the tests read the same three files,
``docs/references.yaml``, ``docs/molbio_facts.yaml`` and ``configs/phix174_genes.yaml``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

from genomeasm.molbio import structure as st

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = ROOT / "docs" / "references.yaml"
FACTS = ROOT / "docs" / "molbio_facts.yaml"
PHIX_GENES = ROOT / "configs" / "phix174_genes.yaml"
PDB_1BNA = ROOT / "data" / "raw" / "1BNA.pdb"

CHECK_LABELS = {
    "crossref": "bibliographic record compared with Crossref (2026-11-14)",
    "file": "data file in this project, SHA-256 checked",
    "software": "compared with Biopython by a test",
    "not-checked": "not checked here (TODO-VERIFY)",
    "abstract": "value read in the paper's abstract (PubMed, 2026-11-14)",
    "corroborated": "classic value; agrees with this program's measurement on PDB 1BNA",
    "textbook": "textbook value, not re-read at the source (TODO-VERIFY)",
    "measured": "measured by this program from the data in data/raw",
}
TOPIC_ORDER = (
    "structure", "genetic material", "replication", "gene expression", "genomes",
    "mutation and biotechnology", "sequencing and assembly", "optimisation", "textbook", "data",
)  # fmt: skip


@lru_cache(maxsize=1)
def references() -> dict[str, dict[str, Any]]:
    data = yaml.safe_load(REFERENCES.read_text(encoding="utf-8"))
    return {str(r["key"]): dict(r) for r in data["references"]}


@lru_cache(maxsize=1)
def facts() -> dict[str, dict[str, Any]]:
    data = yaml.safe_load(FACTS.read_text(encoding="utf-8"))
    return {str(f["key"]): dict(f) for f in data["facts"]}


def cite(key: str) -> str:
    """Short citation, for example ``Watson & Crick 1953`` or ``Drew et al. 1981``."""
    ref = references()[key]
    if ref.get("short"):  # data sources are cited by their accession, not by an author
        return str(ref["short"])
    names = [str(a).rsplit(" ", 1)[0] if " " in str(a) and len(str(a)) < 40 else str(a)
             for a in ref["authors"]]  # fmt: skip
    year = ref["year"] or "n.d."
    if len(names) == 1 and not ref["et_al"]:
        who = names[0]
    elif len(names) == 2 and not ref["et_al"]:
        who = f"{names[0]} & {names[1]}"
    else:
        who = f"{names[0]} et al."
    return f"{who} {year}"


def reference_line(ref: dict[str, Any]) -> str:
    """Full reference in one line (authors, year, title, journal, volume, pages, DOI)."""
    authors = ", ".join(str(a) for a in ref["authors"]) + (" et al." if ref["et_al"] else "")
    year = ref["year"] or "n.d."
    place = str(ref["journal"])
    if ref.get("volume"):
        place += f" {ref['volume']}"
    if ref.get("pages"):
        place += f": {str(ref['pages']).replace('-', '–')}"
    line = f"{authors} ({year}). {ref['title']}. {place}."
    if ref.get("doi"):
        line += f" doi:{ref['doi']}"
    return line


def fact_text(key: str) -> str:
    """Statement of a fact followed by its short citation."""
    fact = facts()[key]
    return f"{fact['statement']} ({cite(str(fact['source']))})"


# ---- the genes of phiX174 ------------------------------------------------------------------------


@dataclass(frozen=True)
class Gene:
    name: str
    parts: tuple[tuple[int, int], ...]  # 1-based, inclusive; two parts when it crosses the origin
    product: str
    protein_id: str

    @property
    def length(self) -> int:
        return sum(end - start + 1 for start, end in self.parts)

    def sequence(self, genome: str) -> str:
        return "".join(genome[start - 1 : end] for start, end in self.parts)

    def frame(self) -> int:
        """Reading frame 0, 1 or 2 of the linear record (position of the first codon modulo 3)."""
        return (self.parts[0][0] - 1) % 3


@lru_cache(maxsize=1)
def phix_annotation() -> dict[str, Any]:
    return dict(yaml.safe_load(PHIX_GENES.read_text(encoding="utf-8")))


def phix_genes() -> list[Gene]:
    return [
        Gene(str(g["name"]), tuple((int(a), int(b)) for a, b in g["parts"]),
             str(g["product"]), str(g["protein_id"]))
        for g in phix_annotation()["genes"]
    ]  # fmt: skip


def gene(name: str) -> Gene:
    return next(g for g in phix_genes() if g.name == name)


def overlapping_pairs(genes: list[Gene]) -> list[tuple[str, str, int]]:
    """Pairs of genes that share letters of the genome, with the number of shared positions."""
    covered = {g.name: {p for a, b in g.parts for p in range(a, b + 1)} for g in genes}
    pairs = []
    for i, first in enumerate(genes):
        for second in genes[i + 1 :]:
            shared = len(covered[first.name] & covered[second.name])
            if shared:
                pairs.append((first.name, second.name, shared))
    return pairs


# ---- measured on the experimental structure ------------------------------------------------------


@lru_cache(maxsize=1)
def structure_1bna() -> st.Structure:
    return st.load_structure(str(PDB_1BNA))


@lru_cache(maxsize=1)
def measured_helix() -> dict[str, float]:
    """Helix parameters measured on PDB 1BNA (see ``structure.helix_geometry`` for the method)."""
    return st.helix_geometry(structure_1bna()).summary()


@lru_cache(maxsize=1)
def backbone_template() -> dict[str, list[float]]:
    return st.backbone_template(structure_1bna())
