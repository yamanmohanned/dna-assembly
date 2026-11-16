"""Scenes for the 3D viewer: the experimental structure and a schematic helix of any sequence.

A scene is a plain dictionary (JSON) that the viewer in ``genomeasm/appkit/web`` draws; all the
geometry is decided here, in tested Python. Coordinates are in angstrom, in a right-handed system
with the helix axis along +y, so a right-handed helix is drawn right-handed.

Two kinds of scene:

* ``atomic_scene``: every atom of the experimental structure (PDB 1BNA), with covalent bonds found
  from distances and the Watson-Crick hydrogen bonds measured between the real atoms;
* ``schematic_scene``: one bead per phosphate, sugar and base for any sequence. The beads are
  placed with the average backbone geometry measured on 1BNA and stacked with the canonical rise
  and helical repeat; it is a teaching model, labelled as such wherever it is shown.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from genomeasm.molbio import structure as st
from genomeasm.molbio.sequence import COMPLEMENT, HYDROGEN_BONDS

COMPONENT_CODE = {"phosphate": 0, "sugar": 1, "base": 2}
SCHEMATIC_RISE = 3.4  # angstrom per base pair (Watson and Crick 1953)
SCHEMATIC_BP_PER_TURN = 10.4  # helical repeat in solution (Wang 1979)


def _round(values: Any, digits: int = 3) -> list[float]:
    return [round(float(v), digits) for v in np.asarray(values, dtype=np.float64).ravel()]


def _extent(xyz: np.ndarray, *, flat: bool) -> list[float]:
    """Half-width, half-height and half-depth the camera must show (angstrom, with a margin).

    A helix is turned about its vertical axis by the viewer, so its half-width is its largest
    distance from that axis; a flat scene (a base pair seen from above) is not turned.
    """
    half_y = float(np.max(np.abs(xyz[:, 1]))) + 3.0
    if flat:
        return [round(float(np.max(np.abs(xyz[:, 0]))) + 3.0, 2), round(half_y, 2),
                round(float(np.max(np.abs(xyz[:, 2]))) + 1.0, 2)]  # fmt: skip
    radial = float(np.max(np.hypot(xyz[:, 0], xyz[:, 2]))) + 2.5
    return [round(radial, 2), round(half_y, 2), round(radial, 2)]


def _helix_basis(geometry: st.HelixGeometry) -> np.ndarray:
    """Rows = the new x, y, z axes: y along the helix axis, z from the middle pair's strand B to
    strand A (so the viewer first looks at the helix from the side of that pair)."""
    y = geometry.axis
    middle = geometry.frames[len(geometry.frames) // 2]
    z = middle.long_axis - np.dot(middle.long_axis, y) * y
    z = z / np.linalg.norm(z)
    x = np.cross(y, z)
    return np.array([x, y, z])


def _residue_record(structure: st.Structure, residue: st.Residue, partner: int | None,
                    n_hbonds: int, strand: int) -> dict[str, Any]:  # fmt: skip
    return {
        "base": residue.base,
        "chain": residue.chain,
        "num": residue.resseq,
        "strand": strand,
        "partner": partner,
        "hb": n_hbonds,
        "name": st.NUCLEOTIDE_NAMES[residue.base],
    }


def atomic_scene(
    structure: st.Structure, residues: list[int] | None = None, orient: str = "helix"
) -> dict[str, Any]:
    """The experimental structure, or only some of its residues (a nucleotide, a base pair).

    ``orient`` is ``"helix"`` (axis vertical) or ``"pair"`` (looking down the helix axis onto the
    plane of the bases, the purine on the left).
    """
    geometry = st.helix_geometry(structure)
    basis = _helix_basis(geometry)
    keep = set(range(len(structure.residues))) if residues is None else set(residues)
    if orient == "pair":
        frame = next(
            (f for f in geometry.frames if f.res_a in keep or f.res_b in keep), geometry.frames[0]
        )
        normal = geometry.axis
        long_axis = frame.long_axis - np.dot(frame.long_axis, normal) * normal
        long_axis = long_axis / np.linalg.norm(long_axis)
        if not structure.residues[frame.res_a].is_purine:
            long_axis = -long_axis  # purine on the left: the line points from pyrimidine to purine
        x = -long_axis
        basis = np.array([x, np.cross(normal, x), normal])
    atom_ids = [i for r in structure.residues if r.index in keep for i in r.atoms]
    new_index = {old: new for new, old in enumerate(atom_ids)}
    centre = structure.coords[atom_ids].mean(axis=0)
    xyz = (structure.coords[atom_ids] - centre) @ basis.T

    partner: dict[int, int] = {}
    n_hb: dict[int, int] = {}
    hbonds = []
    for pair in st.base_pairs(structure):
        partner[pair.purine], partner[pair.pyrimidine] = pair.pyrimidine, pair.purine
        n_hb[pair.purine] = n_hb[pair.pyrimidine] = len(pair.hbonds)
        if pair.purine not in keep or pair.pyrimidine not in keep:
            continue
        bases = (structure.residues[pair.purine].base, structure.residues[pair.pyrimidine].base)
        for bond in pair.hbonds:
            hbonds.append(
                {
                    "a": _round(xyz[new_index[bond.atom_a]]),
                    "b": _round(xyz[new_index[bond.atom_b]]),
                    "d": round(bond.distance, 2),
                    "ra": pair.purine,
                    "rb": pair.pyrimidine,
                    "text": f"{bases[0]} {bond.names[0]} ··· {bond.names[1]} {bases[1]}",
                    "role": st.hbond_description(bases, bond.names),
                }
            )
    first_chain = structure.chains()[0]
    records = [
        _residue_record(structure, r, partner.get(r.index), n_hb.get(r.index, 0),
                        0 if r.chain == first_chain else 1)
        for r in structure.residues
    ]  # fmt: skip
    bonds = [
        (new_index[i], new_index[j])
        for i, j in st.covalent_bonds(structure)
        if i in new_index and j in new_index
    ]
    labels = []
    if residues is None:
        for chain in structure.chains():
            members = [r for r in structure.residues if r.chain == chain]
            for residue, atom, text in ((members[0], "O5'", "5′"), (members[-1], "O3'", "3′")):
                index = structure.atom_index(residue, atom)
                if index is not None:
                    p = xyz[new_index[index]]
                    out = np.array([p[0], 0.0, p[2]])
                    out = out / (np.linalg.norm(out) or 1.0)
                    labels.append({"xyz": _round(p + 2.6 * out + [0, np.sign(p[1]) * 1.2, 0]),
                                   "text": text, "kind": "end",
                                   "strand": 0 if chain == first_chain else 1})  # fmt: skip
        middle = geometry.frames[len(geometry.frames) // 2]
        y_local = middle.long_axis - np.dot(middle.long_axis, geometry.axis) * geometry.axis
        major = np.cross(y_local / np.linalg.norm(y_local), geometry.axis)  # x of the pair frame
        on_axis = middle.centre - centre
        for sign, text in ((1.0, "major groove"), (-1.0, "minor groove")):
            labels.append({"xyz": _round((on_axis + sign * 14.0 * major) @ basis.T),
                           "text": text, "kind": "groove"})  # fmt: skip
    return {
        "kind": "atomic",
        "pdb": structure.pdb_id,
        "atoms": {
            "xyz": _round(xyz),
            "el": [structure.atoms[i].element for i in atom_ids],
            "name": [structure.atoms[i].name for i in atom_ids],
            "comp": [COMPONENT_CODE[structure.atoms[i].component] for i in atom_ids],
            "res": [next(r.index for r in structure.residues if i in r.atoms) for i in atom_ids],
        },
        "bonds": [k for pair_ in bonds for k in pair_],
        "hbonds": hbonds,
        "residues": records,
        "labels": labels,
        "radius": round(float(np.max(np.linalg.norm(xyz, axis=1))) + 2.0, 2),
        "extent": _extent(xyz, flat=orient == "pair"),
    }


def schematic_scene(
    sequence: str,
    template: dict[str, list[float]],
    *,
    rise: float = SCHEMATIC_RISE,
    bp_per_turn: float = SCHEMATIC_BP_PER_TURN,
    first_position: int = 1,
) -> dict[str, Any]:
    """A double helix of ``sequence`` (strand 1, 5' to 3' upwards) and its complement.

    ``template`` is ``structure.backbone_template``: where the phosphate, the sugar and the
    pairing edge of the base sit relative to the base pair. Each nucleotide is three beads:
    phosphate (P), sugar (S) and base (B); bonds P-S and S-B inside a nucleotide and S-P along
    the strand; 2 or 3 hydrogen bonds between paired bases.
    """
    n = len(sequence)
    if n < 2:
        raise ValueError("a helix needs at least 2 base pairs")
    twist = 2.0 * math.pi / bp_per_turn
    xyz: list[list[float]] = []
    el: list[str] = []
    name: list[str] = []
    comp: list[int] = []
    res: list[int] = []
    bonds: list[int] = []
    hbonds: list[dict[str, Any]] = []
    residues: list[dict[str, Any]] = []
    height0 = -(n - 1) * rise / 2.0

    def place(local: list[float], i: int, flip: bool) -> list[float]:
        """Local pair-frame coordinates (x, y, z: z along the axis) to scene coordinates."""
        x, y, z = local
        if flip:  # strand 2: the same nucleotide turned half a turn about the x axis of the pair
            y, z = -y, -z
        angle = i * twist
        # rotation by +angle about the axis (right-handed); pair (x, y, z) -> scene (x', z', y')
        gx = x * math.cos(angle) - y * math.sin(angle)
        gy = x * math.sin(angle) + y * math.cos(angle)
        return [gx, height0 + i * rise + z, -gy]

    index_of: dict[tuple[int, int, str], int] = {}
    for strand in (0, 1):
        for i in range(n):
            base = sequence[i] if strand == 0 else COMPLEMENT[sequence[i]]
            residue_index = len(residues)
            residues.append(
                {
                    "base": base,
                    "chain": "1" if strand == 0 else "2",
                    "num": first_position + i,
                    "strand": strand,
                    "partner": i + n if strand == 0 else i,
                    "hb": HYDROGEN_BONDS[base],
                    "name": st.NUCLEOTIDE_NAMES[base],
                }
            )
            for key, element, label, code in (
                ("P", "P", "phosphate", 0),
                ("sugar", "S", "deoxyribose", 1),
                ("edge", "B", "base", 2),
            ):
                index_of[(strand, i, key)] = len(el)
                xyz.append(place(template[key], i, strand == 1))
                el.append(element)
                name.append(label)
                comp.append(code)
                res.append(residue_index)
            bonds += [index_of[(strand, i, "P")], index_of[(strand, i, "sugar")]]
            bonds += [index_of[(strand, i, "sugar")], index_of[(strand, i, "edge")]]
    for i in range(n - 1):
        # strand 1 runs 5' to 3' upwards: sugar i is linked to the phosphate of nucleotide i + 1;
        # strand 2 runs the other way: sugar i + 1 is linked to the phosphate of nucleotide i
        bonds += [index_of[(0, i, "sugar")], index_of[(0, i + 1, "P")]]
        bonds += [index_of[(1, i + 1, "sugar")], index_of[(1, i, "P")]]
    for i in range(n):
        a = np.array(xyz[index_of[(0, i, "edge")]])
        b = np.array(xyz[index_of[(1, i, "edge")]])
        count = HYDROGEN_BONDS[sequence[i]]
        major = np.array(place([1.0, 0.0, 0.0], i, False)) - np.array(
            place([0.0, 0.0, 0.0], i, False)
        )
        for k in range(count):
            shift = (k - (count - 1) / 2.0) * 1.15 * major
            pair_text = f"{sequence[i]}·{COMPLEMENT[sequence[i]]}"
            hbonds.append(
                {"a": _round(a + shift), "b": _round(b + shift),
                 "d": round(float(np.linalg.norm(a - b)), 2), "ra": i, "rb": i + n,
                 "text": f"{pair_text} hydrogen bond {k + 1} of {count}", "role": ""}
            )  # fmt: skip
    top, bottom = height0 + (n - 1) * rise, height0
    labels = []
    for strand, i, text, dy in ((0, 0, "5′", -3.2), (0, n - 1, "3′", 3.2),
                                (1, 0, "3′", -3.2), (1, n - 1, "5′", 3.2)):  # fmt: skip
        p = np.array(xyz[index_of[(strand, i, "P")]])
        labels.append({"xyz": _round([p[0] * 1.12, (top if dy > 0 else bottom) + dy, p[2] * 1.12]),
                       "text": text, "kind": "end", "strand": strand})  # fmt: skip
    middle = n // 2
    origin = np.array(place([0.0, 0.0, 0.0], middle, False))
    major = np.array(place([1.0, 0.0, 0.0], middle, False)) - origin
    for sign, text in ((1.0, "major groove"), (-1.0, "minor groove")):
        labels.append({"xyz": _round(origin + sign * 15.0 * major), "text": text, "kind": "groove"})
    coords = np.array(xyz)
    return {
        "kind": "schematic",
        "pdb": "",
        "atoms": {"xyz": _round(coords), "el": el, "name": name, "comp": comp, "res": res},
        "bonds": bonds,
        "hbonds": hbonds,
        "residues": residues,
        "labels": labels,
        "radius": round(float(np.max(np.linalg.norm(coords, axis=1))) + 3.0, 2),
        "extent": _extent(coords, flat=False),
        "rise": rise,
        "bp_per_turn": bp_per_turn,
    }
