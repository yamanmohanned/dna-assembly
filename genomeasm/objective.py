"""The objective (fitness) function and the asymmetric travelling-salesman (ATSP) encoding.

Plain-language story: put the reads in some order and glue each read to the one before it. Every
glue saves as many letters as the two reads overlap. The *fitness* F of an order is the total
number of letters saved: bigger is better. The *length* L of the glued text is what remains:
the sum of all read lengths minus F, so maximising F and minimising L are the same goal.

Precise definitions (master prompt, "Mathematical specification"), reads indexed 0..n-1:

* A solution is a permutation ``pi`` of the n reads.
* ``F(pi) = sum_{k=1}^{n-1} ov(pi_k, pi_{k+1})``  to be maximised.
* ``L(pi) = sum_i len_i - F(pi)``  to be minimised: the length of the merged sequence.
* ATSP cost matrix on n + 1 nodes. Node 0 is a dummy; read ``i`` is node ``i + 1``::

      c(i, j) = len_j - ov(i, j)       for reads i != j   (new letters read j adds after read i)
      c(0, j) = len_j                  (the first read adds all its letters)
      c(i, 0) = 0                      (the last read is followed by nothing)

  The diagonal is set to 0 and is never used by a tour.
* Identity (proved in docs/OBJECTIVE_FUNCTION.md and tested): the Hamiltonian cycle
  ``0 -> pi_1 -> ... -> pi_n -> 0`` has cost exactly ``L(pi)``. Hence a minimum-cost tour is an
  order of maximum fitness, and any TSP solver for asymmetric instances can be used.

All costs are integers, as required by LKH-3.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

IntArray = NDArray[np.int64]
IntMatrix = NDArray[np.int32]


def as_order(order: Sequence[int] | IntArray, n: int) -> IntArray:
    """Return ``order`` as an array after checking that it is a permutation of 0..n-1."""
    arr = np.asarray(order, dtype=np.int64)
    if arr.ndim != 1 or arr.size != n or not np.array_equal(np.sort(arr), np.arange(n)):
        raise ValueError(f"order must be a permutation of 0..{n - 1}")
    return arr


def is_permutation(order: Sequence[int] | IntArray, n: int) -> bool:
    arr = np.asarray(order)
    return bool(arr.ndim == 1 and arr.size == n and np.array_equal(np.sort(arr), np.arange(n)))


def fitness(order: Sequence[int] | IntArray, ov: IntMatrix) -> int:
    """F(pi): total overlap along the order (to be maximised)."""
    arr = as_order(order, ov.shape[0])
    if arr.size < 2:
        return 0
    return int(ov[arr[:-1], arr[1:]].sum(dtype=np.int64))


def merged_length(order: Sequence[int] | IntArray, lengths: IntArray, ov: IntMatrix) -> int:
    """L(pi) = sum of read lengths minus F(pi): length of the merged sequence."""
    return int(lengths.sum()) - fitness(order, ov)


def cost_matrix(ov: IntMatrix, lengths: IntArray) -> IntArray:
    """The (n+1) x (n+1) ATSP cost matrix with the dummy node 0 (see module docstring)."""
    n = ov.shape[0]
    if lengths.shape != (n,):
        raise ValueError("lengths must have one entry per read")
    cost = np.zeros((n + 1, n + 1), dtype=np.int64)
    cost[0, 1:] = lengths
    cost[1:, 1:] = lengths[None, :] - ov
    np.fill_diagonal(cost, 0)
    return cost


def tour_from_order(order: Sequence[int] | IntArray) -> IntArray:
    """The cycle 0 -> pi_1 -> ... -> pi_n (the closing edge back to 0 is implied)."""
    arr = np.asarray(order, dtype=np.int64)
    return np.concatenate([[0], arr + 1])


def order_from_tour(tour: Sequence[int] | IntArray, n: int) -> IntArray:
    """Read order of a Hamiltonian cycle on the n + 1 nodes; rotated to start at the dummy node."""
    arr = np.asarray(tour, dtype=np.int64)
    if arr.size == n + 2 and arr[0] == arr[-1]:  # a tour written with the start repeated at the end
        arr = arr[:-1]
    if arr.size != n + 1 or not np.array_equal(np.sort(arr), np.arange(n + 1)):
        raise ValueError(f"tour must visit each of the {n + 1} nodes exactly once")
    start = int(np.flatnonzero(arr == 0)[0])
    rotated = np.roll(arr, -start)
    return rotated[1:] - 1


def cycle_cost(cost: IntArray, tour: Sequence[int] | IntArray) -> int:
    """Cost of the closed tour (the last node returns to the first)."""
    arr = np.asarray(tour, dtype=np.int64)
    return int(cost[arr, np.roll(arr, -1)].sum(dtype=np.int64))


# --------------------------------------------------------------------------------------
# TSPLIB files
# --------------------------------------------------------------------------------------


def write_tsplib_atsp(path: Path, cost: IntArray, *, name: str, comment: str = "") -> None:
    """Write the cost matrix as a TSPLIB ATSP file (FULL_MATRIX, integer weights, LF endings)."""
    if cost.ndim != 2 or cost.shape[0] != cost.shape[1]:
        raise ValueError("cost must be a square matrix")
    if not np.issubdtype(cost.dtype, np.integer):
        raise ValueError("TSPLIB weights must be integers")
    header = [
        f"NAME: {name}",
        "TYPE: ATSP",
        f"COMMENT: {comment}" if comment else "COMMENT: DNA fragment assembly as ATSP",
        f"DIMENSION: {cost.shape[0]}",
        "EDGE_WEIGHT_TYPE: EXPLICIT",
        "EDGE_WEIGHT_FORMAT: FULL_MATRIX",
        "EDGE_WEIGHT_SECTION",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        handle.write("\n".join(header) + "\n")
        np.savetxt(handle, cost, fmt="%d", delimiter=" ", newline="\n")
        handle.write("EOF\n")


def read_tsplib_atsp(path: Path) -> tuple[dict[str, str], IntArray]:
    """Read a TSPLIB ATSP file written by ``write_tsplib_atsp``: (header fields, cost matrix)."""
    text = Path(path).read_text(encoding="ascii")
    fields: dict[str, str] = {}
    numbers: list[str] = []
    in_matrix = False
    for line in text.splitlines():
        line = line.strip()
        if not line or line == "EOF":
            continue
        if line == "EDGE_WEIGHT_SECTION":
            in_matrix = True
        elif in_matrix:
            numbers.extend(line.split())
        elif ":" in line:
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
    if fields.get("TYPE") != "ATSP" or fields.get("EDGE_WEIGHT_FORMAT") != "FULL_MATRIX":
        raise ValueError("only TYPE: ATSP with EDGE_WEIGHT_FORMAT: FULL_MATRIX is supported")
    dimension = int(fields["DIMENSION"])
    if len(numbers) != dimension * dimension:
        raise ValueError(f"expected {dimension * dimension} weights, found {len(numbers)}")
    matrix = np.array(numbers, dtype=np.int64).reshape(dimension, dimension)
    return fields, matrix
