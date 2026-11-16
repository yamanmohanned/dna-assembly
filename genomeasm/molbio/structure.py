"""Atomic structure of DNA: a PDB reader and measurements made on the real coordinates.

The project ships one experimental structure, PDB entry 1BNA (the Drew-Dickerson dodecamer
``CGCGAATTCGCG``, X-ray diffraction at 1.9 angstrom; Drew et al. 1981). Nothing about the helix is
typed in by hand here: base pairs, hydrogen bonds, rise, twist and radius are *measured* from the
atom positions, so the numbers shown in the application are traceable to the deposited file.

Method of the helix measurement (simple on purpose, and stated wherever the numbers are shown):

* the centre of a base pair is the midpoint of the purine C8 and the pyrimidine C6 atom (in B-DNA
  this point lies close to the helix axis);
* the axis is the straight line that fits the centres best (first principal direction);
* rise = distance between neighbouring centres along the axis; twist = angle between the C8/C6
  lines of neighbouring pairs, seen along the axis; a positive twist is a right-handed helix.

A crystal structure is slightly bent and irregular, so a straight axis is an approximation; the
standard deviations are reported with the means.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]

COVALENT_RADIUS = {"C": 0.76, "N": 0.71, "O": 0.66, "P": 1.07, "H": 0.31}  # angstrom
VDW_RADIUS = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "H": 1.20}  # Bondi (1964), angstrom
PHOSPHATE_ATOMS = frozenset({"P", "OP1", "OP2", "O5'", "O3'"})
SUGAR_ATOMS = frozenset({"C1'", "C2'", "C3'", "C4'", "C5'", "O4'"})
BASE_OF_RESIDUE = {"DA": "A", "DC": "C", "DG": "G", "DT": "T"}
NUCLEOTIDE_NAMES = {
    "A": "deoxyadenosine monophosphate (dAMP)",
    "C": "deoxycytidine monophosphate (dCMP)",
    "G": "deoxyguanosine monophosphate (dGMP)",
    "T": "deoxythymidine monophosphate (dTMP)",
}
# Watson-Crick hydrogen bonds as (atom of the purine, atom of the pyrimidine)
WATSON_CRICK = {
    ("A", "T"): (("N6", "O4"), ("N1", "N3")),
    ("G", "C"): (("O6", "N4"), ("N1", "N3"), ("N2", "O2")),
}
HBOND_ROLE = {  # (purine atom, pyrimidine atom) -> which atom carries the hydrogen
    ("N6", "O4"): "N6-H of adenine donates to O4 of thymine",
    ("N1", "N3"): "",  # filled per pair below
    ("O6", "N4"): "N4-H of cytosine donates to O6 of guanine",
    ("N2", "O2"): "N2-H of guanine donates to O2 of cytosine",
}
MAX_HBOND = 3.5  # angstrom between the two heavy atoms; typical N-H...O/N bonds are 2.7 to 3.1


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    xyz: tuple[float, float, float]
    element: str

    @property
    def component(self) -> str:
        """phosphate, sugar or base: the three parts of a nucleotide."""
        if self.name in PHOSPHATE_ATOMS:
            return "phosphate"
        return "sugar" if self.name in SUGAR_ATOMS else "base"


@dataclass(frozen=True)
class Residue:
    index: int  # position in ``Structure.residues``
    chain: str
    resseq: int
    resname: str
    atoms: tuple[int, ...]  # indices into ``Structure.atoms``

    @property
    def base(self) -> str:
        return BASE_OF_RESIDUE[self.resname]

    @property
    def is_purine(self) -> bool:
        return self.base in "AG"


@dataclass(frozen=True)
class HydrogenBond:
    atom_a: int  # atom index (purine side)
    atom_b: int  # atom index (pyrimidine side)
    names: tuple[str, str]
    distance: float


@dataclass(frozen=True)
class BasePair:
    purine: int  # residue index
    pyrimidine: int  # residue index
    hbonds: tuple[HydrogenBond, ...]


@dataclass
class Structure:
    pdb_id: str
    title: str
    experiment: str
    resolution: float | None
    doi: str
    atoms: list[Atom]
    residues: list[Residue]
    n_water: int
    _coords: FloatArray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._coords = np.array([a.xyz for a in self.atoms], dtype=np.float64)

    @property
    def coords(self) -> FloatArray:
        return self._coords

    def chains(self) -> list[str]:
        return sorted({r.chain for r in self.residues})

    def sequence(self, chain: str) -> str:
        return "".join(r.base for r in self.residues if r.chain == chain)

    def atom_index(self, residue: Residue, name: str) -> int | None:
        for i in residue.atoms:
            if self.atoms[i].name == name:
                return i
        return None

    def position(self, residue: Residue, name: str) -> FloatArray:
        index = self.atom_index(residue, name)
        if index is None:
            raise KeyError(f"residue {residue.resname}{residue.resseq} has no atom {name}")
        return self._coords[index]


def parse_pdb(text: str) -> Structure:
    """Read the DNA atoms of a PDB file (fixed columns of the PDB format, version 3)."""
    title, experiment, doi, pdb_id = [], "", "", ""
    resolution: float | None = None
    atoms: list[Atom] = []
    n_water = 0
    for line in text.splitlines():
        record = line[:6].strip()
        if record == "HEADER":
            pdb_id = line[62:66].strip()
        elif record == "TITLE":
            title.append(line[10:80].strip())
        elif record == "EXPDTA":
            experiment = line[10:80].strip()
        elif record == "JRNL" and line[12:16].strip() == "DOI":
            doi = line[19:80].strip()
        elif record == "REMARK" and line[6:10].strip() == "2" and "RESOLUTION." in line:
            resolution = float(line.split("RESOLUTION.")[1].split()[0])
        elif record == "HETATM" and line[17:20].strip() == "HOH":
            n_water += 1
        elif record == "ATOM" and line[17:20].strip() in BASE_OF_RESIDUE:
            if line[16] not in " A":  # alternate locations: keep the first
                continue
            name = line[12:16].strip()
            element = line[76:78].strip() or name[0]
            atoms.append(
                Atom(
                    serial=int(line[6:11]),
                    name=name,
                    resname=line[17:20].strip(),
                    chain=line[21],
                    resseq=int(line[22:26]),
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                    element=element.capitalize(),
                )
            )
    if not atoms:
        raise ValueError("no DNA atoms found: not a PDB file of a DNA structure")
    residues: list[Residue] = []
    start = 0
    for i in range(1, len(atoms) + 1):
        if i == len(atoms) or (atoms[i].chain, atoms[i].resseq) != (
            atoms[start].chain,
            atoms[start].resseq,
        ):
            first = atoms[start]
            residues.append(
                Residue(len(residues), first.chain, first.resseq, first.resname,
                        tuple(range(start, i)))
            )  # fmt: skip
            start = i
    return Structure(
        pdb_id, " ".join(title), experiment, resolution, doi.lower(), atoms, residues, n_water
    )


@lru_cache(maxsize=4)
def load_structure(path: str | Path) -> Structure:
    return parse_pdb(Path(path).read_text(encoding="ascii"))


# ---- chemistry -----------------------------------------------------------------------------------


def covalent_bonds(structure: Structure, tolerance: float = 0.4) -> list[tuple[int, int]]:
    """Bonds found from distances: two atoms closer than the sum of their covalent radii plus a
    tolerance are bonded. Inside a residue every pair is tested; between neighbours of a chain
    only the phosphodiester link O3'(i) to P(i+1)."""
    coords, atoms = structure.coords, structure.atoms
    bonds = []
    for residue in structure.residues:
        ids = residue.atoms
        for a_pos, i in enumerate(ids):
            for j in ids[a_pos + 1 :]:
                limit = (
                    COVALENT_RADIUS[atoms[i].element]
                    + COVALENT_RADIUS[atoms[j].element]
                    + tolerance
                )
                if float(np.linalg.norm(coords[i] - coords[j])) < limit:
                    bonds.append((i, j))
    for first, second in zip(structure.residues, structure.residues[1:], strict=False):
        if first.chain != second.chain:
            continue
        o3, p = structure.atom_index(first, "O3'"), structure.atom_index(second, "P")
        limit = COVALENT_RADIUS["O"] + COVALENT_RADIUS["P"] + tolerance
        if (
            o3 is not None
            and p is not None
            and float(np.linalg.norm(coords[o3] - coords[p])) < limit
        ):
            bonds.append((o3, p))
    return bonds


def base_pairs(structure: Structure) -> list[BasePair]:
    """Watson-Crick pairs found from the coordinates (purine N1 near pyrimidine N3)."""
    coords = structure.coords
    purines = [r for r in structure.residues if r.is_purine]
    pyrimidines = [r for r in structure.residues if not r.is_purine]
    pairs = []
    for purine in purines:
        n1 = structure.position(purine, "N1")
        best: tuple[float, Residue] | None = None
        for other in pyrimidines:
            if (purine.base, other.base) not in WATSON_CRICK or other.chain == purine.chain:
                continue
            distance = float(np.linalg.norm(n1 - structure.position(other, "N3")))
            if distance < MAX_HBOND and (best is None or distance < best[0]):
                best = (distance, other)
        if best is None:
            continue
        other = best[1]
        bonds = []
        for name_a, name_b in WATSON_CRICK[(purine.base, other.base)]:
            a, b = structure.atom_index(purine, name_a), structure.atom_index(other, name_b)
            if a is None or b is None:
                continue
            distance = float(np.linalg.norm(coords[a] - coords[b]))
            if distance < MAX_HBOND:
                bonds.append(HydrogenBond(a, b, (name_a, name_b), distance))
        pairs.append(BasePair(purine.index, other.index, tuple(bonds)))
    return pairs


def hbond_description(pair: tuple[str, str], names: tuple[str, str]) -> str:
    """Which group donates the hydrogen in a Watson-Crick hydrogen bond."""
    if names == ("N1", "N3"):
        if pair == ("A", "T"):
            return "N3-H of thymine donates to N1 of adenine"
        return "N1-H of guanine donates to N3 of cytosine"
    return HBOND_ROLE[names]


# ---- helix geometry ------------------------------------------------------------------------------


@dataclass(frozen=True)
class PairFrame:
    """A base pair seen from the helix: its centre and the direction from strand B to strand A."""

    pair: BasePair
    res_a: int  # residue of the reference strand (the first chain)
    res_b: int
    centre: FloatArray
    long_axis: FloatArray  # unit vector from the strand-B base to the strand-A base, C6/C8 line


@dataclass(frozen=True)
class HelixGeometry:
    n_pairs: int
    origin: FloatArray
    axis: FloatArray  # unit vector, pointing 5' to 3' along the first chain
    frames: tuple[PairFrame, ...]
    rises: tuple[float, ...]
    twists: tuple[float, ...]
    p_radii: tuple[float, ...]
    c1_distances: tuple[float, ...]
    hbond_lengths: tuple[float, ...]

    def summary(self) -> dict[str, float]:
        rise, twist = float(np.mean(self.rises)), float(np.mean(self.twists))
        return {
            "n_pairs": self.n_pairs,
            "rise_mean": rise,
            "rise_sd": float(np.std(self.rises)),
            "twist_mean": twist,
            "twist_sd": float(np.std(self.twists)),
            "bp_per_turn": 360.0 / twist,
            "pitch": rise * 360.0 / twist,
            "p_radius_mean": float(np.mean(self.p_radii)),
            "p_diameter": 2.0 * float(np.mean(self.p_radii)),
            "c1_distance_mean": float(np.mean(self.c1_distances)),
            "hbond_mean": float(np.mean(self.hbond_lengths)),
            "hbond_min": float(np.min(self.hbond_lengths)),
            "hbond_max": float(np.max(self.hbond_lengths)),
            "n_hbonds": len(self.hbond_lengths),
        }


def _ring_atom(structure: Structure, residue: Residue) -> FloatArray:
    return structure.position(residue, "C8" if residue.is_purine else "C6")


def helix_geometry(structure: Structure) -> HelixGeometry:
    """Measure the double helix (method in the module docstring)."""
    pairs = base_pairs(structure)
    first_chain = structure.chains()[0]
    frames = []
    for pair in pairs:
        r1, r2 = structure.residues[pair.purine], structure.residues[pair.pyrimidine]
        res_a, res_b = (r1, r2) if r1.chain == first_chain else (r2, r1)
        pa, pb = _ring_atom(structure, res_a), _ring_atom(structure, res_b)
        direction = pa - pb
        frames.append(
            PairFrame(pair, res_a.index, res_b.index, (pa + pb) / 2.0,
                      direction / np.linalg.norm(direction))
        )  # fmt: skip
    frames.sort(key=lambda f: structure.residues[f.res_a].resseq)
    centres = np.array([f.centre for f in frames])
    origin = centres.mean(axis=0)
    _, _, vt = np.linalg.svd(centres - origin)
    axis = vt[0]
    if float(np.dot(centres[-1] - centres[0], axis)) < 0:
        axis = -axis  # point 5' to 3' along the first chain
    rises, twists = [], []
    for f0, f1 in zip(frames, frames[1:], strict=False):
        rises.append(float(np.dot(f1.centre - f0.centre, axis)))
        v0 = f0.long_axis - np.dot(f0.long_axis, axis) * axis
        v1 = f1.long_axis - np.dot(f1.long_axis, axis) * axis
        angle = np.degrees(np.arctan2(np.dot(np.cross(v0, v1), axis), np.dot(v0, v1)))
        twists.append(float(angle))
    p_radii = []
    for atom, xyz in zip(structure.atoms, structure.coords, strict=True):
        if atom.name == "P":
            offset = xyz - origin
            p_radii.append(float(np.linalg.norm(offset - np.dot(offset, axis) * axis)))
    c1 = [
        float(
            np.linalg.norm(
                structure.position(structure.residues[f.res_a], "C1'")
                - structure.position(structure.residues[f.res_b], "C1'")
            )
        )
        for f in frames
    ]
    lengths = [h.distance for f in frames for h in f.pair.hbonds]
    return HelixGeometry(
        len(frames), origin, axis, tuple(frames), tuple(rises), tuple(twists),
        tuple(p_radii), tuple(c1), tuple(lengths),
    )  # fmt: skip


def backbone_template(structure: Structure) -> dict[str, list[float]]:
    """Average position of the backbone groups of a nucleotide in the frame of its base pair.

    Frame of a pair: z along the helix axis (5' to 3' of the nucleotide's own strand), y from the
    partner base to the nucleotide's base, x = y cross z, origin on the axis at the height of the
    pair. Returned: mean coordinates (angstrom) of the phosphorus atom, the centre of the sugar
    ring, the C1' atom, and the atom of the base that forms the central hydrogen bond (N1 of a
    purine, N3 of a pyrimidine). Used to place the beads of the schematic helix.
    """
    geometry = helix_geometry(structure)
    axis, origin = geometry.axis, geometry.origin
    sums: dict[str, list[FloatArray]] = {"P": [], "sugar": [], "C1'": [], "edge": []}
    for frame in geometry.frames:
        height = float(np.dot(frame.centre - origin, axis))
        on_axis = origin + height * axis
        for residue_index, sign in ((frame.res_a, 1.0), (frame.res_b, -1.0)):
            residue = structure.residues[residue_index]
            z = sign * axis
            y = sign * (frame.long_axis - np.dot(frame.long_axis, axis) * axis)
            y = y / np.linalg.norm(y)
            x = np.cross(y, z)
            basis = np.array([x, y, z])

            def local(point: FloatArray, basis: FloatArray = basis, o: FloatArray = on_axis) -> Any:
                return basis @ (point - o)

            if structure.atom_index(residue, "P") is not None:
                sums["P"].append(local(structure.position(residue, "P")))
            ring = [structure.position(residue, n) for n in ("C1'", "C2'", "C3'", "C4'", "O4'")]
            sums["sugar"].append(local(np.mean(ring, axis=0)))
            sums["C1'"].append(local(structure.position(residue, "C1'")))
            sums["edge"].append(
                local(structure.position(residue, "N1" if residue.is_purine else "N3"))
            )
    return {key: [float(v) for v in np.mean(values, axis=0)] for key, values in sums.items()}
