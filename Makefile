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

.PHONY: help venv test lint format typecheck check app

help:
	@echo Targets: venv test lint format typecheck check
	@echo venv      create .venv and install requirements.lock
	@echo test      run the full pytest suite
	@echo lint      ruff check and ruff format check
	@echo format    apply ruff formatting
	@echo typecheck mypy
	@echo check     lint, typecheck and test
	@echo app       start the interactive application (offline)

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

# TODO (Phase 8): all, which must reproduce every figure, table, report and deck from raw data.

app:
	$(PY) -m genomeasm app run
