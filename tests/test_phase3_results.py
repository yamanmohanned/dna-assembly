"""The committed Phase 3 results: overlap validation, figures and the walking-skeleton runs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from genomeasm import plots, provenance

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "configs" / "phase3_overlap.yaml"
RESULTS_DIR = ROOT / "results" / "phase3"
RESULTS_PATH = RESULTS_DIR / "overlap_validation.json"
FIGURES = ROOT / "reports" / "figures"
STEMS = ["fig02_overlap_length_distribution", "fig03_overlap_heatmap"]


@pytest.fixture(scope="module")
def committed() -> dict[str, Any]:
    return json.loads(RESULTS_PATH.read_text(encoding="utf-8"))


def test_every_check_passed(committed: dict[str, Any]) -> None:
    assert committed["all_passed"] is True and committed["n_failed"] == 0
    assert committed["n_checks"] == len(committed["checks"]) >= 90
    assert all(c["passed"] for c in committed["checks"])


def test_results_match_the_current_config_and_have_full_provenance(
    committed: dict[str, Any],
) -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    assert committed["config"] == cfg, "configs/phase3_overlap.yaml changed: rerun the validation"
    prov = committed["provenance"]
    assert prov["config_hash"] == provenance.config_hash(cfg)
    assert prov["seed"] == cfg["seed"] and len(prov["git"]["commit"]) == 40
    assert prov["git"]["dirty"] is False  # generated at a clean commit
    assert prov["library_versions"]["numba"] and prov["machine"]["cpu_count"]


def test_every_required_validation_is_present(committed: dict[str, Any]) -> None:
    ids = {c["id"] for c in committed["checks"]}
    for prefix in (
        "overlap.recall_error_free",
        "objective.identities",
        "assembly.true_order_contigs_in_genome",
        "overlap.fast_equals_naive",
        "overlap.tolerant_with_zero_tolerance_equals_exact",
        "overlap.preprocessed_reads_are_substring_free",
        "noise.exact_overlaps_lose_true_pairs",
        "noise.tolerant_overlaps_recover_true_pairs",
    ):
        assert any(i.startswith(prefix) for i in ids), f"no check for {prefix}"


def test_error_free_reads_lose_no_true_overlap_at_any_setting(committed: dict[str, Any]) -> None:
    rows = committed["instances"]
    assert len(rows) == 2 * 3 * 7  # genomes x coverages x tau settings
    assert all(row["recall"]["mean"] == 1.0 and row["recall"]["se"] == 0.0 for row in rows)


def test_recommended_min_overlap_gives_an_almost_pure_graph(committed: dict[str, Any]) -> None:
    auto = [row for row in committed["instances"] if row["tau_setting"] == "auto"]
    assert len(auto) == 6
    for row in auto:
        assert row["exact_share_of_edges"]["mean"] >= 0.999
        assert row["spurious"]["mean"] <= 5 * max(1.0, row["spurious_bound"]["mean"])


def test_a_too_small_min_overlap_lets_chance_matches_in(committed: dict[str, Any]) -> None:
    for accession in ("NC_001422.1", "NC_001416.1"):
        rows = {
            (r["coverage"], r["tau_setting"]): r
            for r in committed["instances"]
            if r["accession"] == accession
        }
        assert rows[(10.0, 8)]["spurious"]["mean"] > 10 * max(
            rows[(10.0, 16)]["spurious"]["mean"], 0.1
        )


def test_error_tolerant_overlaps_matter_when_reads_have_errors(committed: dict[str, Any]) -> None:
    modes = committed["noise"]["modes"]
    assert modes["tolerant"]["recall"]["mean"] > modes["exact"]["recall"]["mean"] + 0.5


@pytest.mark.parametrize("stem", STEMS)
def test_figures_exist_with_caption_and_300_dpi(stem: str, committed: dict[str, Any]) -> None:
    from PIL import Image

    for ext in ("svg", "pdf", "png"):
        assert (FIGURES / f"{stem}.{ext}").stat().st_size > 1000
    assert b"\r" not in (FIGURES / f"{stem}.svg").read_bytes()
    sidecar = json.loads((FIGURES / f"{stem}.caption.json").read_text(encoding="utf-8"))
    assert sidecar["config_hash"] == committed["provenance"]["config_hash"]
    assert sidecar["source"] == "results/phase3/overlap_validation.json" and sidecar["caption"]
    with Image.open(FIGURES / f"{stem}.png") as image:
        assert round(image.info["dpi"][0]) == 300


def test_figures_are_redrawn_identically_from_the_saved_results(
    committed: dict[str, Any], tmp_path: Path
) -> None:
    plots.make_phase3_figures(committed, tmp_path, "results/phase3/overlap_validation.json")
    for stem in STEMS:
        assert (tmp_path / f"{stem}.svg").read_bytes() == (FIGURES / f"{stem}.svg").read_bytes()


# ---- walking skeleton on the real phiX174 genome ---------------------------------------------


def _skeleton(name: str) -> dict[str, Any]:
    return json.loads((RESULTS_DIR / name / "summary.json").read_text(encoding="utf-8"))


def test_skeleton_with_the_true_order_recovers_the_genome() -> None:
    summary = _skeleton("skeleton_true")
    assert summary["accession"] == "NC_001422.1" and summary["order"] == "true"
    assert all(summary["contigs_exact_in_genome"])
    assert summary["n_contigs"] <= 2
    assert summary["merged_length_L"] == sum(summary["contig_lengths"])
    assert summary["provenance"]["git"]["dirty"] is False


def test_skeleton_with_the_trivial_order_is_far_worse() -> None:
    true, trivial = _skeleton("skeleton_true"), _skeleton("skeleton_trivial")
    assert true["instance"] == trivial["instance"]  # same reads, same overlap graph
    assert trivial["fitness_F"] * 50 < true["fitness_F"]
    assert trivial["merged_length_L"] > 5 * true["merged_length_L"]
    assert trivial["n_contigs"] > 100
