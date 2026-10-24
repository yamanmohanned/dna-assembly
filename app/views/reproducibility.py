"""Page 8, Reproducibility: configuration, seeds, library versions, commit hash, checksums."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import repro, results_view  # noqa: E402

st.title("Reproducibility: everything needed to repeat the results")
config = repro.config_info()
env = repro.environment()

st.subheader("Code and machine")
c1, c2, c3 = st.columns(3)
commit = env["git"]["commit"] or "unknown"
c1.metric("Git commit of this checkout", commit[:10])
c2.metric(
    "Uncommitted code changes",
    "unknown" if env["git"]["dirty"] is None else ("yes" if env["git"]["dirty"] else "no"),
)
c3.metric("LKH-3 on this computer", "available" if env["lkh3_binary"] else "not available")
st.caption(
    f"Full commit: {commit}. Python {env['machine']['python']} on {env['machine']['platform']}, {env['machine']['cpu_count']} logical CPUs."
)

st.subheader("Configuration and seeds")
st.markdown(
    f"All settings of the experiments are in `{config['path']}` (file SHA-256 `{config['file_sha256'][:16]}...`, "
    f"canonical configuration hash `{config['config_hash'][:16]}...`)."
)
ui.table(
    [
        {"item": "master seed", "value": config["master_seed"]},
        {"item": "read sets per instance setting", "value": config["seeds"]["evaluation"]},
        {"item": "tuning seeds per tuning instance", "value": config["seeds"]["tuning"]},
        {"item": "scaling seeds", "value": config["seeds"]["scaling"]},
        {"item": "time budget per solver (s)", "value": config["budget_seconds"]},
        {"item": "tuning trials per algorithm", "value": config["tuning_trials"]},
    ]
)
st.caption(config["seed_scheme"])

st.subheader("Library versions")
ui.table([{"library": k, "version": v or "not installed"} for k, v in env["libraries"].items()])

st.subheader("Data checksums")
rows = repro.verify_checksums()
ui.table(
    [
        {"file": r["file"], "matches SHA256SUMS": "yes" if r["ok"] else "NO", "sha256": r["actual"]}
        for r in rows
    ]
)

if results_view.available():
    p5 = results_view.load_phase5()
    st.subheader("Where the precomputed results came from")
    st.markdown(
        "Raw runs do not carry a commit, so every stage is recorded in `results/phase5/stages.json` "
        "(*recorded* by the program, or *inferred* from the git history where marked)."
    )
    ui.table(pd.DataFrame(p5.stages["stages"])[["stage", "commit", "dirty", "evidence"]])
    st.markdown("**Analysis provenance** (statistics, tables and figures 4 to 11)")
    st.json(p5.info["provenance"], expanded=False)
    st.subheader("Result files")
    ui.table(repro.result_files())

st.subheader("Repeat it")
st.code(
    "make venv                     # exact packages from requirements.lock\n"
    "make check                    # lint, types and all tests\n"
    "genomeasm experiments tune    # then: run, scale, analyze (resumable)\n"
    "genomeasm app run             # this application",
    language="bash",
)
first, last, count = repro.decision_range()
st.caption(
    f"Every design decision is in docs/DECISIONS.md ({first} to {last}, {count} decisions); "
    "each phase has a report in docs/phases/."
)
