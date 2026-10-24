"""Locations of the project files the application reads (all relative to the repository root)."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
CONFIG = ROOT / "configs" / "phase5_experiments.yaml"
GENOMES_CONFIG = ROOT / "configs" / "genomes.yaml"
GLOSSARY = ROOT / "docs" / "glossary.yaml"
RESULTS = ROOT / "results" / "phase5"
FIGURES = ROOT / "reports" / "figures"
DOCS = ROOT / "docs"

GENOMES = {
    "NC_001422.1": "phiX174 (5,386 letters, a small virus of bacteria)",
    "NC_001416.1": "lambda (48,502 letters, a virus of bacteria)",
}
