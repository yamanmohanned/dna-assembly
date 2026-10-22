# The objective (fitness) function, the solution representation and the ATSP encoding

For the oral exam. Code: [genomeasm/objective.py](../genomeasm/objective.py) (objective, cost matrix,
tours, TSPLIB), [genomeasm/overlap.py](../genomeasm/overlap.py) (overlaps),
[genomeasm/assembly.py](../genomeasm/assembly.py) (layout, consensus),
[genomeasm/instance.py](../genomeasm/instance.py) (everything together). Tests:
`tests/test_overlap.py`, `tests/test_objective.py`, `tests/test_assembly.py`.
Every number in the examples below is produced by the code (toy example) or by
`results/phase3/overlap_validation.json` (real data).

## 1. The story in plain words

A sequencing machine cuts a genome into many short pieces called **reads**, and we must put the
pieces back together. Two reads can be glued when the **end of one equals the beginning of the
other**, like two jigsaw pieces whose edges match. The longer the matching edge, the more letters we
save when we glue them, and the more likely they really were neighbours.

Put the reads in some **order** and glue every read to the one before it. Add up the letters saved
at each glue. That total is the **fitness F** of the order: **bigger is better**. The text we end up
with is shorter by exactly F letters compared with writing all reads one after the other, so
"maximise the savings" and "make the glued text as short as possible" are the same goal.

## 2. The precise definitions

* Reads `r_1, ..., r_n` with lengths `l_1, ..., l_n`, made substring-free (no duplicates, no read
  inside another read).
* **Overlap** `ov(i, j)`: the largest `k` with `tau <= k < min(l_i, l_j)` such that the last `k`
  letters of `r_i` equal the first `k` letters of `r_j`; if there is none, `ov(i, j) = 0`.
  It is **directional**: `ov(i, j)` and `ov(j, i)` are different numbers in general. `tau` is the
  minimum overlap (a parameter, chosen in Section 7).
* **Solution** (the representation): a **permutation** `pi = (pi_1, ..., pi_n)` of the reads. The
  search space has `n!` elements.
* **Fitness**: `F(pi) = sum_{k=1}^{n-1} ov(pi_k, pi_{k+1})`, to be **maximised**.
* **Length of the merged sequence**: `L(pi) = (sum_i l_i) - F(pi)`, to be **minimised**.

### Worked example (the toy example of the master prompt)

Reads A = `ATGCG`, B = `GCGTA`, C = `GTACC`, `tau = 1`. The overlap matrix (rows: first read, columns:
second read) is

|   | A | B | C |
|---|---|---|---|
| **A** | - | 3 | 1 |
| **B** | 1 | - | 3 |
| **C** | 0 | 0 | - |

A then B: the last 3 letters of A (`GCG`) equal the first 3 of B, so gluing saves 3. B then A saves
only 1 (`A`), so the matrix is not symmetric. All six orders (total length of the reads: 15):

| order | glued overlaps | F | L | merged string(s) |
|---|---|---|---|---|
| ABC | 3 + 3 | 6 | 9 | ATGCGTACC |
| BCA | 3 + 0 | 3 | 12 | GCGTACC + ATGCG |
| CAB | 0 + 3 | 3 | 12 | GTACC + ATGCGTA |
| BAC | 1 + 1 | 2 | 13 | GCGTATGCGTACC |
| ACB | 1 + 0 | 1 | 14 | ATGCGTACC + GCGTA |
| CBA | 0 + 1 | 1 | 14 | GTACC + GCGTATGCG |

The best order is A, B, C with `F = 6` and `L = 15 - 6 = 9`: the merged string `ATGCGTACC` has exactly
9 letters. Where two neighbours have overlap 0 the chain breaks and the result is two separate
pieces (contigs), whose lengths still add up to `L`.

## 3. Why `L = sum of lengths - F` (proof)

Glue the reads in the order `pi`. After `pi_1` the text has `l_{pi_1}` letters. When read `pi_{k+1}`
is glued to `pi_k`, its first `ov(pi_k, pi_{k+1})` letters already sit at the end of `pi_k` (this is
the definition of the overlap, and because `ov < l_{pi_k}` the overlap lies inside the last read), so only
`l_{pi_{k+1}} - ov(pi_k, pi_{k+1})` new letters are appended. By induction the text after the last read has

`l_{pi_1} + sum_{k=1}^{n-1} (l_{pi_{k+1}} - ov(pi_k, pi_{k+1})) = sum_i l_i - F(pi)`

letters. Test: `test_merged_length_equals_the_length_of_the_merged_string` builds the merged string
letter by letter on random instances and compares its length with `sum(l) - F`.

## 4. The travelling-salesman encoding and its proof

The **travelling salesman problem (TSP)** asks for the cheapest round trip through a set of towns.
Here the towns are the reads, the trip is the order of the reads, and the price of going from read
`i` to read `j` is **the number of new letters that `j` adds**. Prices depend on the direction
(A to B is not B to A), so this is the **asymmetric TSP (ATSP)**.

Add one extra town, a **dummy node 0**, so that the open chain of reads becomes a closed round trip.
Read `i` is node `i + 1`. The cost matrix `c` on `n + 1` nodes is

* `c(i, j) = l_j - ov(i, j)` for reads `i != j`,
* `c(0, j) = l_j` (the first read adds all its letters),
* `c(i, 0) = 0` (nothing follows the last read),
* diagonal entries are 0 and never used.

All entries are integers (LKH-3 needs integers).

**Claim.** The cycle `0 -> pi_1 -> ... -> pi_n -> 0` costs exactly `L(pi)`.

**Proof.** The cost is `c(0, pi_1) + sum_{k=1}^{n-1} c(pi_k, pi_{k+1}) + c(pi_n, 0)
= l_{pi_1} + sum_{k=1}^{n-1} (l_{pi_{k+1}} - ov(pi_k, pi_{k+1})) + 0 = sum_i l_i - F(pi) = L(pi)`. A round trip that
passes through 0 and every read exactly once, read in the direction of travel starting at 0, lists the reads in a
unique order, so round trips and permutations correspond one to one. Hence **a cheapest round trip is
an order of maximum fitness** (and a shortest merged sequence). Any ATSP solver can therefore be used. QED

Tests: `test_cycle_cost_with_the_dummy_node_equals_merged_length` (random real-like instances),
`test_the_identity_holds_for_any_overlap_matrix_and_lengths` (the proof uses nothing about the
letters, so it holds for any matrix; tested on random matrices), a brute-force comparison over all
orders for small `n`, and the identity on real data in all 42 settings of the Phase 3 validation
(35 with the full cost matrix, the rest through layout and plain merge because the dense
matrix would not fit in memory).

For the toy example the matrix (rows and columns: dummy, A, B, C) is

```
0 5 5 5
0 0 2 4
0 4 0 2
0 5 5 0
```

and the tour `0 -> A -> B -> C -> 0` costs `5 + 2 + 2 + 0 = 9 = L`.

## 5. Why this is hard

There are `n!` orders. For the 516 reads of phiX174 at 10x coverage that is a number with more than
1,000 digits. The ATSP is NP-hard in general (a standard fact; the citation is `TODO-VERIFY`), so
the project compares heuristics (greedy, tabu search, ALNS) with a strong reference solver (LKH-3).
Exact methods work only for tiny `n`: brute force up to 9 reads, Held-Karp dynamic programming up to
15 (Phase 4).

## 6. What the objective does not know

* **Repeats.** If the genome contains the same stretch twice, an order that jumps between the two
  copies can have a higher `F` than the true order. The optimum of `F` is therefore not guaranteed to
  be the true genome order. The project reports both `F` and agreement with the ground truth, and
  `F` of the true ordering as a reference line (declared simplification).
* **Errors.** With the error-tolerant overlap, mismatches are allowed but not charged: the score of an
  overlap is its length. Insertions and deletions are not handled by the tolerant overlap (Hamming
  distance only).
* **Real contigs.** An order always chains all reads. Edges with overlap 0 are the places where the
  assembly is cut into contigs.

## 7. Choosing the minimum overlap `tau`

Too small a `tau` lets unrelated reads match by chance; too large a `tau` throws away true
overlaps and splits contigs. For two unrelated random reads the chance of a match of at least `tau`
letters is below `(4/3) * 4^(-tau)`, so among `n(n-1)` ordered pairs fewer than
`n(n-1) * (4/3) * 4^(-tau)` chance overlaps are expected. `recommend_min_overlap` takes the smallest
`tau` for which this is below 1. Measured on real genomes (Phase 3 report): with that rule, 0.2 to
1.6 chance edges per instance among 1,171 to 137,208 edges, and every true overlap found.

## 8. Where each piece lives

| Idea | Function |
|---|---|
| overlap matrix | `overlap_matrix` (fast), `overlap_matrix_naive` (definition), `overlap_matrix_tolerant` |
| make reads substring-free | `make_substring_free` |
| `F`, `L` | `objective.fitness`, `objective.merged_length` |
| ATSP matrix, tours | `objective.cost_matrix`, `tour_from_order`, `order_from_tour`, `cycle_cost` |
| file for LKH-3 | `objective.write_tsplib_atsp` (`AssemblyInstance.write_tsplib`) |
| gluing and voting | `assembly.layout`, `assembly.consensus`, `assembly.assemble` |

Delta evaluation (updating `F` in constant time after a swap or insert move) is added with the
solvers in Phase 4 and tested there against a full recomputation.
