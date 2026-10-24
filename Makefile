# GNU make. Recipes are single plain commands so they work with sh or cmd.
# Not yet runnable on this machine until GNU make is installed (see docs/phases/PHASE_0.md).
# Every target is a thin wrapper around a python -m command that can also be typed by hand.

ifeq ($(OS),Windows_NT)
PYTHON ?= python
PY := .venv/Scripts/python.exe
else
PYTHON ?= python3
PY := .venv/bin/python
endif

.PHONY: help venv test lint format typecheck check app docs verify-data reproduce all experiments

help:
	@echo Targets: venv test lint format typecheck check
	@echo venv      create .venv and install requirements.lock
	@echo test      run the full pytest suite
	@echo lint      ruff check and ruff format check
	@echo format    apply ruff formatting
	@echo typecheck mypy
	@echo check     lint, typecheck and test
	@echo app       start the interactive application (offline)
	@echo docs      build the report, the slides and the oral-exam sheet
	@echo all       venv if missing, rebuild every table, figure and document from results/, then check
	@echo experiments  rerun the experiments of Phase 5 (hours; not part of all)

venv:
	$(PYTHON) -m venv .venv
	$(PY) -m pip install -r requirements.lock
	$(PY) -m pip install -e . --no-deps

test:
	$(PY) -m pytest

lint:
	$(PY) -m ruff check .
	$(PY) -m ruff format --check .

format:
	$(PY) -m ruff format .

# Phase 1: package and tests. The app folder is added when it has code (Phase 6).
typecheck:
	$(PY) -m mypy genomeasm tests

check: lint typecheck test


app:
	$(PY) -m genomeasm app run

docs:
	$(PY) -m genomeasm docs build

# Phase 8. all: from a clean checkout to every table, figure, report, deck and audit.
# The raw runs of Phase 5 (results/phase5/*.jsonl.gz) are committed; rerunning them takes hours.
$(PY):
	$(PYTHON) -m venv .venv
	$(PY) -m pip install -r requirements.lock
	$(PY) -m pip install -e . --no-deps

verify-data:
	$(PY) -m genomeasm data verify

reproduce: verify-data
	$(PY) -m genomeasm simulate figures
	$(PY) -m genomeasm overlap figures
	$(PY) -m genomeasm experiments analyze
	$(PY) -m genomeasm docs build
	$(PY) -m genomeasm docs verify

all: $(PY) reproduce check

experiments:
	$(PY) -m genomeasm experiments tune
	$(PY) -m genomeasm experiments run
	$(PY) -m genomeasm experiments scale
	$(PY) -m genomeasm experiments analyze
