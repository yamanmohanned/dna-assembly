"""Layout, consensus, instances from simulated reads, and the walking-skeleton command."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from genomeasm import cli
from genomeasm.assembly import assemble, consensus, layout, merged_string
from genomeasm.instance import AssemblyInstance
from genomeasm.io import Genome, NcbiRecord, download_genome, load_genome, parse_fasta
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads

ROOT = Path(__file__).resolve().parent.parent
TOY = ["ATGCG", "GCGTA", "GTACC"]


def toy() -> AssemblyInstance:
    return AssemblyInstance.from_sequences(TOY, min_overlap=1)


# ---- layout -------------------------------------------------------------------------------------


def test_toy_layout_offsets() -> None:
    inst = toy()
    groups, offsets = layout([0, 1, 2], inst.lengths, inst.ov)
    assert groups == [[0, 1, 2]] and offsets == [[0, 2, 4]]  # A at 0, B at 5-3, C at 2+5-3


def test_chain_breaks_where_the_overlap_is_zero() -> None:
    inst = toy()
    groups, offsets = layout([0, 2, 1], inst.lengths, inst.ov)  # ov(A,C)=1, ov(C,B)=0
    assert groups == [[0, 2], [1]] and offsets == [[0, 4], [0]]
    result = inst.assemble([0, 2, 1])
    assert result.n_contigs == 2
    assert result.sequences == ["ATGCGTACC", "GCGTA"]
    assert result.total_length == inst.merged_length([0, 2, 1]) == 14  # 15 letters minus F = 1


def test_plain_merge_of_the_toy_example() -> None:
    inst = toy()
    assert merged_string([0, 1, 2], inst.reads, inst.ov) == "ATGCGTACC"
    # B then A saves 1 letter, A then C saves 1 letter: GCGTA + TGCG + TACC
    assert merged_string([1, 0, 2], inst.reads, inst.ov) == "GCGTATGCGTACC"
    assert merged_string([], [], np.zeros((0, 0), dtype=np.int32)) == ""


# ---- consensus ----------------------------------------------------------------------------------


def test_consensus_outvotes_a_single_error() -> None:
    reads = ["ACGTACGT", "ACGAACGT", "ACGTACGT"]  # the middle read has an error at column 3
    assert consensus(reads, [0, 1, 2], [0, 0, 0]) == "ACGTACGT"


def test_consensus_over_staggered_reads() -> None:
    reads = ["ACGTAC", "GTACGG", "ACGGTT"]
    assert consensus(reads, [0, 1, 2], [0, 2, 4]) == "ACGTACGGTT"


def test_consensus_ties_go_to_the_first_letter_of_acgt() -> None:
    assert consensus(["T", "C"], [0, 1], [0, 0]) == "C"
    assert consensus(["G", "A"], [0, 1], [0, 0]) == "A"


def test_contig_lengths_always_add_up_to_the_merged_length() -> None:
    rng = np.random.default_rng(3)
    genome = "".join(rng.choice(list("ACGT"), size=400))
    reads = [genome[s : s + 40] for s in rng.integers(0, 360, size=60)]
    inst = AssemblyInstance.from_sequences(reads, min_overlap=6)
    for _ in range(10):
        order = rng.permutation(inst.n)
        assert inst.assemble(order).total_length == inst.merged_length(order)


# ---- instances built from simulated reads -------------------------------------------------------


@pytest.fixture(scope="module")
def phix_reads() -> tuple[Genome, AssemblyInstance]:
    genome = load_genome("NC_001422.1", ROOT / "data" / "raw").genome
    reads = simulate_reads(
        genome, SimulationConfig(coverage=10, read_length=ReadLengthModel.fixed(100), seed=21)
    )
    return genome, AssemblyInstance.from_readset(reads)


def test_instance_summary_and_bookkeeping(phix_reads: tuple[Genome, AssemblyInstance]) -> None:
    _, inst = phix_reads
    summary = inst.summary()
    assert summary["n_input_reads"] == 539
    assert (
        summary["n_input_reads"] - summary["n_duplicates_removed"] - summary["n_contained_removed"]
        == summary["n_reads"]
        == inst.n
    )
    assert inst.ov.shape == (inst.n, inst.n) and np.all(np.diag(inst.ov) == 0)
    assert inst.min_overlap == 10  # recommended for about 500 reads
    assert inst.true_starts is not None and inst.true_starts.shape == (inst.n,)


def test_gluing_real_reads_in_true_order_reproduces_the_genome(
    phix_reads: tuple[Genome, AssemblyInstance],
) -> None:
    genome, inst = phix_reads
    order = inst.true_order()
    result = inst.assemble(order)
    assert result.n_contigs <= 2
    assert all(contig.sequence in genome.sequence for contig in result.contigs)
    # the true order is far better than the order in which the reads came out of the sequencer
    assert inst.fitness(order) > 100 * inst.fitness(np.arange(inst.n))
    assert inst.merged_length(order) < 1.01 * genome.length


def test_instance_without_ground_truth_has_no_true_order() -> None:
    with pytest.raises(ValueError, match="no ground truth"):
        toy().true_order()


def test_consensus_repairs_substitution_errors_of_real_reads() -> None:
    genome = load_genome("NC_001422.1", ROOT / "data" / "raw").genome
    reads = simulate_reads(
        genome,
        SimulationConfig(
            coverage=12,
            read_length=ReadLengthModel.fixed(100),
            substitution_rate=0.01,
            seed=44,
        ),
    )
    inst = AssemblyInstance.from_readset(reads, min_overlap=20, tolerance=0.05)
    assert inst.true_starts is not None
    result = assemble(inst.true_order(), inst.reads, inst.lengths, inst.ov)
    wrong = total = wrong_deep = total_deep = 0
    for contig in result.contigs:
        start = int(inst.true_starts[contig.reads[0]])
        window = genome.sequence[start : start + len(contig.sequence)]
        assert len(window) == len(contig.sequence)
        depth = np.zeros(len(window), dtype=np.int64)
        for member, offset in zip(contig.reads, contig.offsets, strict=True):
            depth[offset : offset + len(inst.reads[member])] += 1
        differs = np.array([a != b for a, b in zip(contig.sequence, window, strict=True)])
        deep = depth >= 5
        wrong, total = wrong + int(differs.sum()), total + differs.size
        wrong_deep, total_deep = wrong_deep + int(differs[deep].sum()), total_deep + int(deep.sum())
    raw_rate = float(np.sum(reads.n_substitutions)) / float(np.sum(reads.spans))
    assert raw_rate > 0.009  # the input really is noisy (1% substitutions)
    # Where at least 5 reads vote, a wrong consensus letter needs 3 wrong votes for the same
    # letter: probability about 1e-6 per column, so at most 1 wrong column is allowed here.
    assert total_deep > 3000 and wrong_deep <= 1
    # Tolerant overlaps miss some true neighbours, so the chain breaks into many short contigs
    # whose edge columns have few votes. Voting is still never worse than the raw reads.
    assert wrong / total < raw_rate


# ---- the walking skeleton (command line) --------------------------------------------------------


@pytest.fixture
def fake_cache(tmp_path: Path, ncbi_record: Callable[..., NcbiRecord]) -> Path:
    raw = tmp_path / "raw"
    record = ncbi_record("NC_999999.1")
    download_genome("NC_999999.1", raw, fetcher=lambda a: record)
    return raw


@pytest.mark.parametrize("order", ["trivial", "random", "true"])
def test_cli_assemble_runs_end_to_end(
    order: str, fake_cache: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / f"out_{order}"
    code = cli.main(
        [
            "assemble", "NC_999999.1", "--coverage", "6", "--read-length", "40",
            "--order", order, "--seed", "3", "--raw-dir", str(fake_cache), "--out-dir", str(out),
        ]
    )  # fmt: skip
    assert code == 0
    summary = json.loads((out / "summary.json").read_text())
    assert summary["order"] == order
    assert summary["merged_length_L"] == sum(summary["contig_lengths"])
    assert summary["instance"]["n_reads"] <= summary["instance"]["n_input_reads"]
    assert summary["provenance"]["seed"] == 3 and len(summary["provenance"]["config_hash"]) == 64
    contigs = parse_fasta((out / "contigs.fasta").read_text())
    assert [len(c.sequence) for c in contigs] == summary["contig_lengths"]
    if order == "true":
        assert all(summary["contigs_exact_in_genome"])
    assert "contig(s)" in capsys.readouterr().out


def test_cli_true_order_beats_trivial_order(fake_cache: Path, tmp_path: Path) -> None:
    results = {}
    for order in ("trivial", "true"):
        out = tmp_path / order
        cli.main(
            ["assemble", "NC_999999.1", "--coverage", "8", "--read-length", "40", "--order", order,
             "--seed", "5", "--raw-dir", str(fake_cache), "--out-dir", str(out)]
        )  # fmt: skip
        results[order] = json.loads((out / "summary.json").read_text())
    assert results["true"]["fitness_F"] > results["trivial"]["fitness_F"]
    assert results["true"]["merged_length_L"] < results["trivial"]["merged_length_L"]


def test_cli_assemble_reports_a_missing_genome(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = cli.main(["assemble", "NC_999999.1", "--coverage", "5", "--raw-dir", str(tmp_path)])
    assert code == 1 and "error:" in capsys.readouterr().err
