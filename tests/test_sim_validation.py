"""The validation machinery itself, the committed Phase 2 results and the figure files."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml
from hypothesis import given
from hypothesis import strategies as st

from genomeasm import coverage as cov
from genomeasm import plots, provenance
from genomeasm.io import load_genome
from genomeasm.sim_validation import edit_distance, histogram_z_scores, run_validation
from genomeasm.simulate import ReadLengthModel, SimulationConfig, sample_geometry

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "configs" / "phase2_validation.yaml"
RESULTS_PATH = ROOT / "results" / "phase2" / "validation.json"
FIGURES = ROOT / "reports" / "figures"


# ---- edit distance -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [("", "", 0), ("abc", "", 3), ("", "abc", 3), ("kitten", "sitting", 3), ("ACGT", "ACGT", 0)],
)
def test_edit_distance_known_values(a: str, b: str, expected: int) -> None:
    assert edit_distance(a, b) == expected


@given(a=st.text(alphabet="ACGT", max_size=12), b=st.text(alphabet="ACGT", max_size=12))
def test_edit_distance_is_a_symmetric_metric_bounded_by_lengths(a: str, b: str) -> None:
    d = edit_distance(a, b)
    assert d == edit_distance(b, a)
    assert abs(len(a) - len(b)) <= d <= max(len(a), len(b))
    assert (d == 0) == (a == b)


# ---- the histogram check has teeth ---------------------------------------------------------------


def _histograms(genome_length: int, coverage: float, runs: int, seed: int) -> tuple[Any, ...]:
    cfg = SimulationConfig(coverage=coverage, read_length=ReadLengthModel.fixed(100))
    n_reads = cfg.n_reads(genome_length)
    support, pmf = cfg.read_length.pmf()
    bins = 40
    hists = np.empty((runs, bins))
    for r, child in enumerate(np.random.SeedSequence(seed).spawn(runs)):
        starts, spans = sample_geometry(genome_length, cfg, np.random.default_rng(child))
        hists[r] = cov.depth_histogram(cov.coverage_depth(genome_length, starts, spans), bins)
    return hists, n_reads, support, pmf, bins


def test_histogram_check_accepts_the_true_model_and_rejects_a_biased_one() -> None:
    g = 5386
    hists, n_reads, support, pmf, bins = _histograms(g, 5, runs=200, seed=17)
    truth = cov.expected_depth_histogram(g, n_reads, support, pmf, bins)
    _, z_true = histogram_z_scores(hists, truth, g, 100)
    assert np.max(np.abs(z_true)) < 5.0
    # negative controls: pretend 10% more reads, or a shifted mean depth
    more = cov.expected_depth_histogram(g, round(n_reads * 1.1), support, pmf, bins)
    _, z_more = histogram_z_scores(hists, more, g, 100)
    assert np.max(np.abs(z_more)) > 5.0
    fewer = cov.expected_depth_histogram(g, round(n_reads * 0.9), support, pmf, bins)
    _, z_fewer = histogram_z_scores(hists, fewer, g, 100)
    assert np.max(np.abs(z_fewer)) > 5.0


def test_histogram_check_handles_zero_spread_groups() -> None:
    hists = np.tile(np.array([0.5, 0.3, 0.2]), (10, 1))
    groups, z = histogram_z_scores(hists, np.array([0.5, 0.3, 0.2]), 1000, 100)
    assert groups and np.all(z == 0)
    _, z_bad = histogram_z_scores(hists, np.array([0.4, 0.4, 0.2]), 1000, 100)
    assert np.isinf(z_bad).any()


# ---- the whole pipeline on a reduced configuration (smoke test; evidence is results/) ------------


def test_reduced_validation_run_passes_every_check() -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    cfg.update(coverages=[2, 5], histogram_coverages=[5], replicates=40)
    cfg["large_genome_check"].update(genome_length=300_000, coverages=[1, 2], replicates=30)
    # The z-variance range [0.6, 1.4] is only meaningful with about 100 runs (the sample variance
    # of 15 z-scores varies by +-0.4), so the error checks keep their full number of runs.
    cfg["length_check"].update(genome_length=300_000, n_reads=60_000)
    cfg["substitution_rates"] = [0.01]
    genomes = {
        acc: load_genome(acc, ROOT / "data" / "raw").genome
        for acc in dict.fromkeys([*cfg["genomes"], cfg["error_check"]["genome"]])
    }
    results = run_validation(cfg, genomes, repo_dir=ROOT)
    failed = [c for c in results["checks"] if not c["passed"]]
    assert not failed, failed
    assert results["n_checks"] == len(results["checks"]) > 30
    assert results["provenance"]["config_hash"] == provenance.config_hash(cfg)


# ---- committed results and figures ---------------------------------------------------------------


@pytest.fixture(scope="module")
def committed() -> dict[str, Any]:
    return json.loads(RESULTS_PATH.read_text(encoding="utf-8"))


def test_committed_results_all_passed(committed: dict[str, Any]) -> None:
    assert committed["all_passed"] is True and committed["n_failed"] == 0
    assert committed["n_checks"] == len(committed["checks"]) >= 90
    assert all(c["passed"] for c in committed["checks"])


def test_committed_results_match_the_current_config(committed: dict[str, Any]) -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    assert committed["config"] == cfg, (
        "configs/phase2_validation.yaml changed: rerun the validation"
    )
    assert committed["provenance"]["config_hash"] == provenance.config_hash(cfg)


def test_committed_results_carry_full_provenance(committed: dict[str, Any]) -> None:
    prov = committed["provenance"]
    assert prov["seed"] == committed["config"]["seed"]
    assert len(prov["git"]["commit"]) == 40
    assert prov["library_versions"]["numpy"] and prov["machine"]["cpu_count"]
    assert prov["timestamp_utc"].endswith("+00:00")


def test_committed_results_cover_every_required_validation(committed: dict[str, Any]) -> None:
    ids = {c["id"] for c in committed["checks"]}
    for prefix in (
        "coverage.identity",
        "coverage.uncovered_vs_exact",
        "coverage.histogram_vs_exact",
        "coverage.uncovered_vs_lw_large_genome",
        "contigs.vs_exact",
        "errors.substitution_rate",
        "errors.substitution_calibration",
        "errors.hamming_equals_recorded",
        "errors.deletion_rate",
        "errors.insertion_rate",
        "errors.read_length_identity",
        "errors.edit_distance_bound",
        "lengths.distribution",
        "starts.uniform",
    ):
        assert any(i.startswith(prefix) for i in ids), f"no check for {prefix}"


@pytest.mark.parametrize("stem", ["fig01a_coverage_histogram", "figS1_validation_curves"])
def test_figures_exist_in_all_formats_with_a_caption(stem: str, committed: dict[str, Any]) -> None:
    for ext in ("svg", "pdf", "png"):
        assert (FIGURES / f"{stem}.{ext}").stat().st_size > 1000
    sidecar = json.loads((FIGURES / f"{stem}.caption.json").read_text(encoding="utf-8"))
    assert sidecar["config_hash"] == committed["provenance"]["config_hash"]
    assert sidecar["source"] == "results/phase2/validation.json"
    assert sidecar["caption"]


@pytest.mark.parametrize("stem", ["fig01a_coverage_histogram", "figS1_validation_curves"])
def test_svg_figures_use_unix_line_endings(stem: str) -> None:
    assert b"\r" not in (FIGURES / f"{stem}.svg").read_bytes()


@pytest.mark.parametrize("stem", ["fig01a_coverage_histogram", "figS1_validation_curves"])
def test_png_figures_are_300_dpi(stem: str) -> None:
    from PIL import Image

    with Image.open(FIGURES / f"{stem}.png") as image:
        dpi = image.info["dpi"]
        assert round(dpi[0]) == 300 and round(dpi[1]) == 300
        assert image.width >= 3000


def test_figures_are_redrawn_identically_from_the_saved_results(
    committed: dict[str, Any], tmp_path: Path
) -> None:
    plots.make_phase2_figures(committed, tmp_path, "results/phase2/validation.json")
    for name in ("fig01a_coverage_histogram.svg", "figS1_validation_curves.svg"):
        assert (tmp_path / name).read_bytes() == (FIGURES / name).read_bytes()


# ---- provenance ----------------------------------------------------------------------------------


def test_config_hash_ignores_key_order_but_not_values() -> None:
    assert provenance.config_hash({"a": 1, "b": [1, 2]}) == provenance.config_hash(
        {"b": [1, 2], "a": 1}
    )
    assert provenance.config_hash({"a": 1}) != provenance.config_hash({"a": 2})
    assert len(provenance.config_hash({})) == 64


def test_collect_has_every_required_field() -> None:
    prov = provenance.collect({"x": 1}, seed=5, repo_dir=ROOT)
    assert set(prov) == {
        "config_hash",
        "git",
        "seed",
        "library_versions",
        "machine",
        "timestamp_utc",
    }
    assert prov["seed"] == 5 and prov["library_versions"]["numpy"]


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_git_info_ignores_generated_outputs_but_not_code(tmp_path: Path) -> None:
    def git(*args: str) -> None:
        subprocess.run(
            ["git", "-c", "user.name=t", "-c", "user.email=t@example.org", *args],
            cwd=tmp_path,
            check=True,
            capture_output=True,
        )

    git("init", "-q")
    (tmp_path / "code.py").write_text("x = 1\n")
    git("add", "-A")
    git("commit", "-q", "-m", "init")
    assert provenance.git_info(tmp_path)["dirty"] is False
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "out.json").write_text("{}")
    assert provenance.git_info(tmp_path)["dirty"] is False
    (tmp_path / "code.py").write_text("x = 2\n")
    info = provenance.git_info(tmp_path)
    assert info["dirty"] is True and len(info["commit"]) == 40


def test_git_info_outside_a_repository_is_unknown(tmp_path: Path) -> None:
    info = provenance.git_info(tmp_path)
    assert info == {"commit": None, "dirty": None}
