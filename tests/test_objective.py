"""The objective function and the ATSP encoding: identities (b) and (c) of the master prompt."""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genomeasm import objective
from genomeasm.assembly import layout, merged_string
from genomeasm.instance import AssemblyInstance
from genomeasm.overlap import overlap_matrix

TOY = ["ATGCG", "GCGTA", "GTACC"]


def random_instance(
    rng: np.random.Generator, n_reads: int = 9, genome_length: int = 90
) -> AssemblyInstance:
    """Substring-free reads cut from a random synthetic genome (code-invariant tests only)."""
    genome = "".join(rng.choice(list("ACGT"), size=genome_length))
    reads = []
    for _ in range(n_reads * 3):
        length = int(rng.integers(8, 26))
        start = int(rng.integers(0, genome_length - length + 1))
        reads.append(genome[start : start + length])
    inst = AssemblyInstance.from_sequences(reads, min_overlap=3)
    return inst


# ---- (b) L = sum of lengths - F equals the length of the really merged string -------------------


@settings(deadline=None, max_examples=120)
@given(seed=st.integers(min_value=0, max_value=10**6))
def test_merged_length_equals_the_length_of_the_merged_string(seed: int) -> None:
    rng = np.random.default_rng(seed)
    inst = random_instance(rng)
    assert inst.n >= 2
    for _ in range(4):
        order = rng.permutation(inst.n)
        merged = merged_string(order, inst.reads, inst.ov)
        assert len(merged) == inst.merged_length(order)
        assert len(merged) == int(inst.lengths.sum()) - inst.fitness(order)


@settings(deadline=None, max_examples=80)
@given(seed=st.integers(min_value=0, max_value=10**6))
def test_every_read_sits_in_the_merged_string_at_its_layout_offset(seed: int) -> None:
    rng = np.random.default_rng(seed)
    inst = random_instance(rng)
    order = rng.permutation(inst.n)
    groups, offsets = layout(order, inst.lengths, inst.ov)
    assembly = inst.assemble(order)
    for contig, members, offs in zip(assembly.contigs, groups, offsets, strict=True):
        for member, offset in zip(members, offs, strict=True):
            read = inst.reads[member]
            assert contig.sequence[offset : offset + len(read)] == read  # exact overlaps agree


def test_fitness_of_empty_and_single_orders_is_zero() -> None:
    assert objective.fitness([], np.zeros((0, 0), dtype=np.int32)) == 0
    assert objective.fitness([0], np.zeros((1, 1), dtype=np.int32)) == 0


# ---- (c) the Hamiltonian cycle through the dummy node costs exactly L ---------------------------


@settings(deadline=None, max_examples=150)
@given(seed=st.integers(min_value=0, max_value=10**6))
def test_cycle_cost_with_the_dummy_node_equals_merged_length(seed: int) -> None:
    rng = np.random.default_rng(seed)
    inst = random_instance(rng)
    cost = inst.cost()
    for _ in range(5):
        order = rng.permutation(inst.n)
        tour = objective.tour_from_order(order)
        assert objective.cycle_cost(cost, tour) == inst.merged_length(order)


@settings(deadline=None, max_examples=150)
@given(
    n=st.integers(min_value=1, max_value=9),
    seed=st.integers(min_value=0, max_value=10**6),
)
def test_the_identity_holds_for_any_overlap_matrix_and_lengths(n: int, seed: int) -> None:
    """The proof in docs/OBJECTIVE_FUNCTION.md uses only the definition of c, not the reads."""
    rng = np.random.default_rng(seed)
    lengths = rng.integers(5, 40, size=n).astype(np.int64)
    ov = rng.integers(0, 5, size=(n, n)).astype(np.int32)
    np.fill_diagonal(ov, 0)
    cost = objective.cost_matrix(ov, lengths)
    order = rng.permutation(n)
    expected = int(lengths.sum()) - objective.fitness(order, ov)
    assert objective.cycle_cost(cost, objective.tour_from_order(order)) == expected


def test_maximising_f_and_minimising_the_tour_cost_pick_the_same_order() -> None:
    rng = np.random.default_rng(11)
    for _ in range(12):
        inst = random_instance(rng, n_reads=3, genome_length=60)
        n = min(inst.n, 6)
        sub = AssemblyInstance.from_sequences(inst.reads[:n], min_overlap=3)
        cost = sub.cost()
        orders = list(itertools.permutations(range(sub.n)))
        by_f = max(sub.fitness(o) for o in orders)
        by_cost = min(objective.cycle_cost(cost, objective.tour_from_order(o)) for o in orders)
        assert by_cost == int(sub.lengths.sum()) - by_f


def test_toy_cost_matrix_structure() -> None:
    inst = AssemblyInstance.from_sequences(TOY, min_overlap=1)
    expected = [[0, 5, 5, 5], [0, 0, 2, 4], [0, 4, 0, 2], [0, 5, 5, 0]]
    assert inst.cost().tolist() == expected
    # row 0: the dummy node pays the full length of the first read; column 0: free return
    assert inst.cost()[0, 1:].tolist() == [5, 5, 5] and inst.cost()[1:, 0].tolist() == [0, 0, 0]


def test_costs_are_positive_integers_for_substring_free_reads() -> None:
    inst = random_instance(np.random.default_rng(2))
    cost = inst.cost()
    assert cost.dtype == np.int64
    assert np.all(cost[0, 1:] >= 1)  # the first read adds all its letters
    assert np.all(cost[1:, 1:][~np.eye(inst.n, dtype=bool)] >= 1)  # ov < read length


def test_cost_matrix_checks_its_input() -> None:
    with pytest.raises(ValueError, match="one entry per read"):
        objective.cost_matrix(np.zeros((3, 3), dtype=np.int32), np.array([5, 5], dtype=np.int64))


# ---- orders, tours ------------------------------------------------------------------------------


@given(st.permutations(list(range(7))))
def test_order_and_tour_are_inverse_representations(order: list[int]) -> None:
    tour = objective.tour_from_order(order)
    assert tour[0] == 0 and sorted(tour.tolist()) == list(range(8))
    assert objective.order_from_tour(tour, 7).tolist() == order


def test_a_rotated_or_closed_tour_gives_the_same_order() -> None:
    order = [2, 0, 3, 1]
    tour = objective.tour_from_order(order)  # [0, 3, 1, 4, 2]
    assert objective.order_from_tour(np.roll(tour, 2), 4).tolist() == order
    assert objective.order_from_tour(np.append(tour, 0), 4).tolist() == order


@pytest.mark.parametrize("tour", [[0, 1, 1, 2], [0, 1, 2], [1, 2, 3, 4], [0, 1, 2, 5]])
def test_invalid_tours_are_rejected(tour: list[int]) -> None:
    with pytest.raises(ValueError, match="exactly once"):
        objective.order_from_tour(tour, 3)


@pytest.mark.parametrize("order", [[0, 1, 1], [0, 1], [0, 1, 3], [[0, 1, 2]]])
def test_invalid_orders_are_rejected(order: list[int]) -> None:
    assert not objective.is_permutation(order, 3)
    with pytest.raises(ValueError, match="permutation"):
        objective.fitness(order, overlap_matrix(TOY, 1))


# ---- TSPLIB files -------------------------------------------------------------------------------


def test_tsplib_roundtrip_and_header(tmp_path: Path) -> None:
    inst = AssemblyInstance.from_sequences(TOY, min_overlap=1)
    path = tmp_path / "toy.atsp"
    inst.write_tsplib(path, "toy")
    fields, matrix = objective.read_tsplib_atsp(path)
    assert np.array_equal(matrix, inst.cost())
    assert fields["TYPE"] == "ATSP" and fields["DIMENSION"] == "4"
    assert (
        fields["EDGE_WEIGHT_TYPE"] == "EXPLICIT" and fields["EDGE_WEIGHT_FORMAT"] == "FULL_MATRIX"
    )
    raw = path.read_bytes()
    assert b"\r" not in raw and raw.endswith(b"EOF\n")
    lines = raw.decode().splitlines()
    assert lines[0] == "NAME: toy" and "EDGE_WEIGHT_SECTION" in lines
    assert lines[lines.index("EDGE_WEIGHT_SECTION") + 1] == "0 5 5 5"


def test_tsplib_rejects_non_integer_or_non_square_matrices(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="integers"):
        objective.write_tsplib_atsp(
            tmp_path / "x.atsp",
            np.zeros((3, 3)),  # type: ignore[arg-type]  # floats on purpose
            name="x",
        )
    with pytest.raises(ValueError, match="square"):
        objective.write_tsplib_atsp(tmp_path / "x.atsp", np.zeros((3, 4), dtype=np.int64), name="x")


def test_reading_an_unsupported_or_damaged_file_fails(tmp_path: Path) -> None:
    other = tmp_path / "other.tsp"
    other.write_text("NAME: x\nTYPE: TSP\nDIMENSION: 2\nEOF\n")
    with pytest.raises(ValueError, match="only TYPE: ATSP"):
        objective.read_tsplib_atsp(other)
    short = tmp_path / "short.atsp"
    short.write_text(
        "NAME: x\nTYPE: ATSP\nDIMENSION: 3\nEDGE_WEIGHT_TYPE: EXPLICIT\n"
        "EDGE_WEIGHT_FORMAT: FULL_MATRIX\nEDGE_WEIGHT_SECTION\n0 1 2\nEOF\n"
    )
    with pytest.raises(ValueError, match="expected 9 weights"):
        objective.read_tsplib_atsp(short)
