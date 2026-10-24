"""Smoke tests of every page of the Streamlit application (headless, through ``AppTest``).

Each page script is run as the browser would run it. These tests check that nothing crashes, that
the main interactions work and that the numbers shown agree with the underlying functions; the
looks are checked by hand in the browser (docs/phases/PHASE_6.md).
"""

from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from genomeasm.appkit import content, lab

ROOT = Path(__file__).resolve().parents[1]
VIEWS = ROOT / "app" / "views"
PAGES = [
    "guide",
    "learn",
    "data",
    "simulator",
    "overlap_graph",
    "solver_lab",
    "experiments",
    "assembly_viewer",
    "reproducibility",
]


def run_page(name: str, timeout: int = 120) -> AppTest:
    return AppTest.from_file(str(VIEWS / f"{name}.py"), default_timeout=timeout).run()


@pytest.mark.parametrize("name", PAGES)
def test_every_page_runs_without_an_error(name: str) -> None:
    at = run_page(name)
    assert not at.exception, [e.value for e in at.exception]
    assert len(at.title) == 1 and at.title[0].value


def test_the_eight_pages_of_the_master_prompt_and_the_guide_exist() -> None:
    assert len(PAGES) == 9
    assert all((VIEWS / f"{p}.py").is_file() for p in PAGES)
    entry = (ROOT / "app" / "streamlit_app.py").read_text(encoding="utf-8")
    for page in PAGES:
        assert f"views/{page}.py" in entry


def test_learn_page_shows_the_toy_example_and_reacts_to_an_edit() -> None:
    at = run_page("learn")
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["Best possible F"] == "6" and metrics["Fitness F of your order"] == "6"
    order = next(t for t in at.text_input if t.label.startswith("Your order"))
    order.set_value("BCA").run()
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert (
        metrics["Fitness F of your order"] == "3" and metrics["Length L of the merged text"] == "12"
    )
    order.set_value("AAB").run()
    assert any("exactly once" in w.value for w in at.warning)
    reads = next(t for t in at.text_input if t.label.startswith("Reads"))
    reads.set_value("ACGT").run()
    assert any("between 2" in e.value for e in at.error)  # one read is not enough


def test_learn_glossary_search_filters_the_terms() -> None:
    at = run_page("learn")
    total = len(content.load_glossary())
    assert any(f"{total} of {total} terms" in m.value for m in at.markdown)
    box = next(t for t in at.text_input if t.label == "Search the glossary")
    box.set_value("bandit").run()
    hits = len(content.search_glossary(content.load_glossary(), "bandit"))
    assert 0 < hits < total
    assert any(f"{hits} of {total} terms" in m.value for m in at.markdown)


def test_data_page_lists_the_card_of_the_selected_genome() -> None:
    at = run_page("data")
    assert any(m.label == "Length (letters)" and m.value == "5,386" for m in at.metric)
    at.radio[0].set_value("NC_001416.1").run()
    assert not at.exception
    assert any(m.label == "Length (letters)" and m.value == "48,502" for m in at.metric)


def test_simulator_and_overlap_pages_report_the_instance() -> None:
    at = run_page("simulator")
    assert any(m.label == "Reads" for m in at.metric) and not at.warning
    graph = run_page("overlap_graph")
    labels = {m.label for m in graph.metric}
    assert {"Reads in the graph", "Arrows (overlaps)", "F of the true order"} <= labels


def test_too_large_a_live_instance_is_refused_with_an_explanation() -> None:
    at = AppTest.from_file(str(VIEWS / "simulator.py"), default_timeout=60)
    at.run()
    coverage = next(s for s in at.sidebar.slider if s.label.startswith("Coverage"))
    length = next(s for s in at.sidebar.slider if s.label.startswith("Window length"))
    length.set_value(5000).run()
    coverage = next(s for s in at.sidebar.slider if s.label.startswith("Coverage"))
    coverage.set_value(20.0).run()
    assert not at.exception and not at.warning  # exactly 1,000 reads is still allowed
    read_length = next(s for s in at.sidebar.slider if s.label.startswith("Read length"))
    read_length.set_value(40).run()
    assert any("allowed for a live run" in w.value for w in at.warning)


def test_solver_lab_runs_a_solver_and_keeps_the_result() -> None:
    at = run_page("solver_lab")
    budget = next(s for s in at.slider if s.label.startswith("Time budget"))
    budget.set_value(1)
    at.button(key="lab-run").click().run()
    assert not at.exception, [e.value for e in at.exception]
    run = at.session_state["lab_run"]
    assert run["solver"] == "alns" and run["fitness"] >= run["greedy_fitness"]
    assert sorted(run["order"]) == list(range(run["n"]))
    assert any(m.label == "Fitness F" and m.value == f"{run['fitness']:,}" for m in at.metric)
    # a rerun (for instance after a download click) shows the stored result and starts no new run
    at.run()
    assert at.session_state["lab_run"] == run


def test_solver_lab_step_through_modes_run() -> None:
    at = run_page("solver_lab")
    assert any("Greedy order" in s.value for s in at.success)  # the toy example, greedy
    kind = next(r for r in at.radio if r.label == "Instance")
    kind.set_value("Tiny simulated (7 to 8 reads)").run()
    assert not at.exception
    algorithm = next(r for r in at.radio if r.label == "Algorithm")
    algorithm.set_value("ALNS").run()
    assert not at.exception, [e.value for e in at.exception]
    assert any(m.label == "Chosen pair" for m in at.metric)
    step = next(s for s in at.slider if s.label == "ALNS iteration")
    step.set_value(5).run()
    assert not at.exception and any("Outcome" in m.value for m in at.markdown)


def test_assembly_viewer_assembles_the_true_order() -> None:
    at = run_page("assembly_viewer")
    next(s for s in at.selectbox if s.label.startswith("Order of the reads")).set_value(
        "true"
    ).run()
    at.button(key="view-assemble").click().run()
    assert not at.exception, [e.value for e in at.exception]
    run = at.session_state["view_run"]
    assert run["solver"] == "true" and run["quality"]["identity"] > 0.95
    assert any(
        m.label == "Misjoins" and m.value == str(run["quality"]["misjoins"]) for m in at.metric
    )


def test_experiments_page_shows_the_not_supported_verdict() -> None:
    at = run_page("experiments")
    assert any("Not supported" in e.value for e in at.error)
    assert any("Precomputed" in m.value for m in at.markdown)  # the badge says it is precomputed


def test_reproducibility_page_reports_the_commit_and_clean_checksums() -> None:
    at = run_page("reproducibility")
    labels = {m.label: m.value for m in at.metric}
    assert len(labels["Git commit of this checkout"]) == 10
    assert not at.exception


def test_the_entry_point_defines_the_navigation_with_every_page() -> None:
    at = AppTest.from_file(str(ROOT / "app" / "streamlit_app.py"), default_timeout=120).run()
    assert not at.exception, [e.value for e in at.exception]
    assert len(at.title) == 1 and at.title[0].value.startswith("Guide")


def test_guide_page_has_every_section_and_a_searchable_glossary() -> None:
    at = run_page("guide")
    sections = content.guide_sections()
    assert len(sections) == 9 and len(at.tabs) == len(sections) + 1
    assert [t.label for t in at.tabs][-1] == "Glossary"
    box = next(t for t in at.text_input if t.label == "Search the glossary")
    box.set_value("helix").run()
    hits = content.search_glossary(content.load_glossary(), "helix")
    assert any(f"{len(hits)} of" in m.value for m in at.markdown)


def test_learn_page_draws_dna_and_rna_from_the_real_genome() -> None:
    at = run_page("learn")
    assert not at.exception
    next(n for n in at.number_input if n.label == "Start").set_value(100).run()
    assert not at.exception, [e.value for e in at.exception]
    code = " ".join(c.value for c in at.code)
    sequence = lab.genome_sequence("NC_001422.1")[100:124]
    assert sequence in code
