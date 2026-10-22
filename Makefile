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

.PHONY: help venv test lint format typecheck check

help:
	@echo Targets: venv test lint format typecheck check
	@echo venv      create .venv and install requirements.lock
	@echo test      run the full pytest suite
	@echo lint      ruff check and ruff format check
	@echo format    apply ruff formatting
	@echo typecheck mypy
	@echo check     lint, typecheck and test

venv:
	$(PYTHON) -m venv .venv
	$(PY) -m pip install -r requirements.lock

test:
	$(PY) -m pytest

lint:
	$(PY) -m ruff check .
	$(PY) -m ruff format --check .

format:
	$(PY) -m ruff format .

# Phase 0: only tests exist. Phase 1 extends this to genomeasm and app.
typecheck:
	$(PY) -m mypy tests

check: lint typecheck test

# TODO (Phase 8): all, which must reproduce every figure, table, report and deck from raw data.
