"""Simulator: ground truth, error model, determinism and input validation.

Real data (rule R1) is used wherever the numbers matter. Tiny synthetic genomes appear only in
property tests that check invariants of the code itself, never results.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genomeasm import cli
from genomeasm.io import Genome, NcbiRecord, download_genome, load_genome, parse_fasta
from genomeasm.sim_validation import edit_distance
from genomeasm.simulate import (
    ReadLengthModel,
    ReadSet,
    SimulationConfig,
    sample_geometry,
    simulate_reads,
    streams,
)

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"


@pytest.fixture(scope="module")
def phix() -> Genome:
    return load_genome("NC_001422.1", RAW).genome


def fixed(length: int = 100) -> ReadLengthModel:
    return ReadLengthModel.fixed(length)


def window(genome: Genome, reads: ReadSet, i: int) -> str:
    return genome.sequence[reads.starts[i] : reads.starts[i] + reads.spans[i]]


def is_subsequence(small: str, big: str) -> bool:
    it = iter(big)
    return all(ch in it for ch in small)


# ---- read-length models -------------------------------------------------------------------


def test_fixed_model() -> None:
    support, pmf = fixed(100).pmf()
    assert support.tolist() == [100] and pmf.tolist() == [1.0]
    assert fixed(100).mean == 100


def test_uniform_model() -> None:
    model = ReadLengthModel.uniform(80, 120)
    support, pmf = model.pmf()
    assert support[0] == 80 and support[-1] == 120 and len(support) == 41
    assert pmf == pytest.approx(np.full(41, 1 / 41))
    assert model.mean == pytest.approx(100.0)


def test_normal_model_is_a_proper_truncated_distribution() -> None:
    model = ReadLengthModel.normal(100, 15, 50, 150)
    support, pmf = model.pmf()
    assert pmf.sum() == pytest.approx(1.0)
    assert model.mean == pytest.approx(100.0, abs=1e-6)  # symmetric around mu
    assert pmf[np.searchsorted(support, 100)] == np.max(pmf)


@pytest.mark.parametrize(
    "build",
    [
        lambda: ReadLengthModel("weird", 1, 2),
        lambda: ReadLengthModel.uniform(0, 5),
        lambda: ReadLengthModel.uniform(10, 5),
        lambda: ReadLengthModel("fixed", 5, 6),
        lambda: ReadLengthModel.normal(100, 0, 50, 150),
    ],
)
def test_invalid_length_models_are_rejected(build: Callable[[], ReadLengthModel]) -> None:
    with pytest.raises(ValueError):
        build()


def test_length_model_dict_roundtrip() -> None:
    for model in (fixed(90), ReadLengthModel.uniform(5, 9), ReadLengthModel.normal(50, 5, 30, 70)):
        assert ReadLengthModel.from_dict(model.to_dict()) == model
    assert ReadLengthModel.from_dict({"length": 100}) == fixed(100)


# ---- configuration --------------------------------------------------------------------------


def test_read_count_and_nominal_coverage() -> None:
    cfg = SimulationConfig(coverage=5, read_length=fixed(100))
    assert cfg.n_reads(5386) == 269
    assert cfg.nominal_coverage(5386) == pytest.approx(269 * 100 / 5386)
    assert SimulationConfig(coverage=0.001, read_length=fixed(100)).n_reads(5386) == 1


@pytest.mark.parametrize(
    "kwargs",
    [
        {"coverage": 0},
        {"coverage": -1},
        {"coverage": 5, "substitution_rate": -0.1},
        {"coverage": 5, "insertion_rate": 0.6},
        {"coverage": 5, "deletion_rate": 1.0},
    ],
)
def test_invalid_configs_are_rejected(kwargs: dict[str, float]) -> None:
    with pytest.raises(ValueError):
        SimulationConfig(read_length=fixed(), **kwargs)  # type: ignore[arg-type]


def test_reads_longer_than_the_genome_are_rejected() -> None:
    with pytest.raises(ValueError, match="longer than the genome"):
        SimulationConfig(coverage=5, read_length=fixed(500)).n_reads(100)


def test_config_dict_roundtrip() -> None:
    cfg = SimulationConfig(
        coverage=7.5,
        read_length=ReadLengthModel.normal(100, 10, 70, 130),
        substitution_rate=0.01,
        insertion_rate=0.002,
        deletion_rate=0.003,
        seed=42,
    )
    assert SimulationConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))) == cfg


# ---- ground truth on the real genome ---------------------------------------------------------


def test_error_free_reads_are_exact_substrings_at_the_recorded_position(phix: Genome) -> None:
    reads = simulate_reads(phix, SimulationConfig(coverage=10, read_length=fixed(), seed=1))
    assert reads.n_reads == 539
    for i in range(reads.n_reads):
        assert reads.sequences[i] == window(phix, reads, i)
    assert reads.n_substitutions.sum() == reads.n_insertions.sum() == reads.n_deletions.sum() == 0
    assert np.min(reads.starts) >= 0
    assert np.max(reads.starts + reads.spans) <= phix.length
    assert reads.realised_coverage == pytest.approx(539 * 100 / phix.length)


def test_variable_length_reads_stay_inside_the_genome(phix: Genome) -> None:
    model = ReadLengthModel.normal(100, 20, 40, 160)
    reads = simulate_reads(phix, SimulationConfig(coverage=8, read_length=model, seed=2))
    assert np.min(reads.spans) >= 40 and np.max(reads.spans) <= 160
    assert len(set(reads.spans.tolist())) > 20
    assert np.max(reads.starts + reads.spans) <= phix.length
    assert all(reads.sequences[i] == window(phix, reads, i) for i in range(reads.n_reads))


def test_ids_and_fasta_carry_the_ground_truth(phix: Genome) -> None:
    cfg = SimulationConfig(coverage=3, read_length=fixed(80), substitution_rate=0.02, seed=5)
    reads = simulate_reads(phix, cfg)
    assert reads.ids[0] == "read_000000" and len(set(reads.ids)) == reads.n_reads
    records = parse_fasta(reads.to_fasta())
    assert [r.sequence for r in records] == reads.sequences
    first = records[0].header.split()
    assert first[0] == "read_000000"
    assert f"start={reads.starts[0]}" in first and f"span={reads.spans[0]}" in first
    assert f"sub={reads.n_substitutions[0]}" in first


# ---- determinism --------------------------------------------------------------------------------


def test_same_seed_gives_identical_reads_and_different_seed_does_not(phix: Genome) -> None:
    cfg = SimulationConfig(
        coverage=5, read_length=fixed(), substitution_rate=0.01, insertion_rate=0.005, seed=11
    )
    a, b = simulate_reads(phix, cfg), simulate_reads(phix, cfg)
    assert a.sequences == b.sequences
    assert np.array_equal(a.starts, b.starts) and np.array_equal(a.n_insertions, b.n_insertions)
    other = SimulationConfig.from_dict({**cfg.to_dict(), "seed": 12})
    assert simulate_reads(phix, other).sequences != a.sequences


def test_changing_only_the_error_rate_does_not_move_the_reads(phix: Genome) -> None:
    clean = simulate_reads(phix, SimulationConfig(coverage=5, read_length=fixed(), seed=9))
    noisy = simulate_reads(
        phix, SimulationConfig(coverage=5, read_length=fixed(), substitution_rate=0.05, seed=9)
    )
    assert np.array_equal(clean.starts, noisy.starts)
    assert np.array_equal(clean.spans, noisy.spans)
    assert clean.sequences != noisy.sequences


def test_geometry_sampler_matches_the_full_simulator(phix: Genome) -> None:
    cfg = SimulationConfig(coverage=4, read_length=ReadLengthModel.uniform(60, 90), seed=3)
    starts, spans = sample_geometry(phix.length, cfg, streams(cfg.seed)[0])
    reads = simulate_reads(phix, cfg)
    assert np.array_equal(starts, reads.starts) and np.array_equal(spans, reads.spans)


# ---- error model ---------------------------------------------------------------------------------


def test_substitutions_change_exactly_the_recorded_bases(phix: Genome) -> None:
    cfg = SimulationConfig(coverage=10, read_length=fixed(), substitution_rate=0.05, seed=4)
    reads = simulate_reads(phix, cfg)
    assert reads.n_substitutions.sum() > 0
    for i in range(reads.n_reads):
        original = window(phix, reads, i)
        hamming = sum(x != y for x, y in zip(reads.sequences[i], original, strict=True))
        assert len(reads.sequences[i]) == reads.spans[i]
        assert hamming == reads.n_substitutions[i]  # a substitution always changes the base


def test_insertions_only_add_bases(phix: Genome) -> None:
    cfg = SimulationConfig(coverage=5, read_length=fixed(), insertion_rate=0.05, seed=6)
    reads = simulate_reads(phix, cfg)
    assert reads.n_insertions.sum() > 0 and reads.n_deletions.sum() == 0
    for i in range(reads.n_reads):
        assert len(reads.sequences[i]) == reads.spans[i] + reads.n_insertions[i]
        assert is_subsequence(window(phix, reads, i), reads.sequences[i])


def test_deletions_only_remove_bases(phix: Genome) -> None:
    cfg = SimulationConfig(coverage=5, read_length=fixed(), deletion_rate=0.05, seed=7)
    reads = simulate_reads(phix, cfg)
    assert reads.n_deletions.sum() > 0 and reads.n_insertions.sum() == 0
    for i in range(reads.n_reads):
        assert len(reads.sequences[i]) == reads.spans[i] - reads.n_deletions[i]
        assert is_subsequence(reads.sequences[i], window(phix, reads, i))


def test_mixed_errors_respect_length_identity_and_edit_bound(phix: Genome) -> None:
    cfg = SimulationConfig(
        coverage=4,
        read_length=fixed(),
        substitution_rate=0.03,
        insertion_rate=0.02,
        deletion_rate=0.02,
        seed=8,
    )
    reads = simulate_reads(phix, cfg)
    assert np.array_equal(reads.lengths, reads.spans - reads.n_deletions + reads.n_insertions)
    for i in range(0, reads.n_reads, 9):
        edits = reads.n_substitutions[i] + reads.n_insertions[i] + reads.n_deletions[i]
        assert edit_distance(reads.sequences[i], window(phix, reads, i)) <= edits


def test_reads_that_lose_every_base_are_reported() -> None:
    genome = Genome.from_sequence("NC_000001.1", "ACGT" * 5)
    cfg = SimulationConfig(coverage=50, read_length=fixed(1), deletion_rate=0.5, seed=1)
    with pytest.raises(ValueError, match="lost all its bases"):
        simulate_reads(genome, cfg)


# ---- property test on tiny synthetic genomes (code invariants only) -----------------------------


@settings(deadline=None, max_examples=60)
@given(
    sequence=st.text(alphabet="ACGT", min_size=60, max_size=240),
    read_length=st.integers(min_value=10, max_value=50),
    coverage=st.floats(min_value=0.5, max_value=6),
    sub=st.floats(min_value=0, max_value=0.3),
    ins=st.floats(min_value=0, max_value=0.2),
    dele=st.floats(min_value=0, max_value=0.2),
    seed=st.integers(min_value=0, max_value=2**31),
)
def test_invariants_hold_for_any_small_genome(
    sequence: str,
    read_length: int,
    coverage: float,
    sub: float,
    ins: float,
    dele: float,
    seed: int,
) -> None:
    genome = Genome.from_sequence("NC_000001.1", sequence)
    cfg = SimulationConfig(
        coverage=coverage,
        read_length=fixed(read_length),
        substitution_rate=sub,
        insertion_rate=ins,
        deletion_rate=dele,
        seed=seed,
    )
    reads = simulate_reads(genome, cfg)
    assert reads.n_reads == cfg.n_reads(genome.length)
    assert np.min(reads.starts) >= 0 and np.max(reads.starts + reads.spans) <= genome.length
    assert np.all(reads.spans == read_length)
    assert set("".join(reads.sequences)) <= set("ACGT")
    assert np.array_equal(reads.lengths, reads.spans - reads.n_deletions + reads.n_insertions)
    again = simulate_reads(genome, cfg)
    assert again.sequences == reads.sequences
    if sub == ins == dele == 0:
        assert all(reads.sequences[i] == window(genome, reads, i) for i in range(reads.n_reads))


# ---- command line --------------------------------------------------------------------------------


def test_cli_simulate_run_writes_reads_and_a_config_record(
    tmp_path: Path, ncbi_record: Callable[..., NcbiRecord], capsys: pytest.CaptureFixture[str]
) -> None:
    raw = tmp_path / "raw"
    record = ncbi_record("NC_999999.1")
    download_genome("NC_999999.1", raw, fetcher=lambda a: record)
    out = tmp_path / "out"
    code = cli.main(
        [
            "simulate", "run", "NC_999999.1", "--coverage", "3", "--read-length", "40",
            "--read-length-sd", "5", "--sub-rate", "0.01", "--seed", "7",
            "--raw-dir", str(raw), "--out-dir", str(out),
        ]
    )  # fmt: skip
    assert code == 0
    summary = json.loads((out / "reads.config.json").read_text())
    assert summary["simulation"]["seed"] == 7
    assert summary["simulation"]["read_length"]["kind"] == "normal"
    assert summary["provenance"]["config_hash"] and summary["provenance"]["seed"] == 7
    assert len(parse_fasta((out / "reads.fasta").read_text())) == summary["n_reads"]
    assert b"\r" not in (out / "reads.fasta").read_bytes()
    assert "reads (coverage" in capsys.readouterr().out


def test_cli_simulate_run_reports_a_missing_genome(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = cli.main(
        ["simulate", "run", "NC_999999.1", "--coverage", "3", "--raw-dir", str(tmp_path)]
    )
    assert code == 1
    assert "error:" in capsys.readouterr().err
