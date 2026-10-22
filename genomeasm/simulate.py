"""Sequencing simulator: cut a real genome into reads, with errors, and remember the truth.

Plain-language story: a sequencing machine cannot read a whole genome at once. It reads many short
pieces (reads) from random places, and sometimes makes typing mistakes. This module plays the
machine. Because we know the genome, every simulated read keeps a note of exactly where it came
from and which mistakes were made. That note is the ground truth used later to judge assemblies.

Model (every parameter is explicit and stored in ``SimulationConfig``; rule R2):

* One strand only, reads forward, never reaching outside the genome (decision D-017).
* Number of reads ``N = round(coverage * G / mean_read_length)``, so the nominal coverage is
  ``c = N * mean_read_length / G``.
* Each read draws its length from the read-length model, then a start position uniformly among
  the ``G - length + 1`` positions where the read fits.
* Errors act on the reference window of a read, base by base and independently:
  - with probability ``deletion_rate`` the base is dropped;
  - otherwise, with probability ``substitution_rate`` it is replaced by one of the three other
    bases (each with probability 1/3);
  - after every reference base (also a deleted one), with probability ``insertion_rate`` one
    random base (A, C, G or T with probability 1/4 each) is inserted.
  So the substitution rate is the rate among bases that were kept, and at most one insertion
  follows each reference base. All rates must lie in [0, 0.5].
* Randomness: ``numpy.random.default_rng`` seeded through ``SeedSequence(seed)``; one child
  stream for read geometry (lengths, starts) and one for errors. The same seed gives identical
  reads, and changing only an error rate never moves the read positions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.stats import norm

from genomeasm.io import Genome

IntArray = NDArray[np.int64]
FloatArray = NDArray[np.float64]

_ENCODE = np.zeros(256, dtype=np.uint8)
for _code, _letter in enumerate("ACGT"):
    _ENCODE[ord(_letter)] = _code
_DECODE = np.frombuffer(b"ACGT", dtype=np.uint8)

MAX_ERROR_RATE = 0.5


@dataclass(frozen=True)
class ReadLengthModel:
    """Distribution of read lengths: fixed, discrete uniform, or truncated discretised normal."""

    kind: str
    low: int
    high: int
    mu: float = 0.0
    sigma: float = 0.0

    def __post_init__(self) -> None:
        if self.kind not in {"fixed", "uniform", "normal"}:
            raise ValueError(f"unknown read-length kind {self.kind!r}")
        if self.low < 1 or self.high < self.low:
            raise ValueError("need 1 <= low <= high")
        if self.kind == "fixed" and self.low != self.high:
            raise ValueError("a fixed length needs low == high")
        if self.kind == "normal" and self.sigma <= 0:
            raise ValueError("a normal model needs sigma > 0")

    @classmethod
    def fixed(cls, length: int) -> ReadLengthModel:
        return cls("fixed", length, length)

    @classmethod
    def uniform(cls, low: int, high: int) -> ReadLengthModel:
        return cls("uniform", low, high)

    @classmethod
    def normal(cls, mu: float, sigma: float, low: int, high: int) -> ReadLengthModel:
        """Normal(mu, sigma) rounded to whole numbers and cut to the range [low, high]."""
        return cls("normal", low, high, mu=mu, sigma=sigma)

    def pmf(self) -> tuple[IntArray, FloatArray]:
        """Possible lengths and their probabilities (they sum to 1)."""
        support = np.arange(self.low, self.high + 1, dtype=np.int64)
        if self.kind == "normal":
            edges = np.arange(self.low, self.high + 2) - 0.5
            cdf = norm.cdf(edges, loc=self.mu, scale=self.sigma)
            probs = np.diff(cdf)
            probs = probs / probs.sum()
        else:
            probs = np.full(support.size, 1.0 / support.size)
        return support, probs.astype(np.float64)

    @property
    def mean(self) -> float:
        """Exact mean of the length distribution."""
        support, probs = self.pmf()
        return float(np.dot(support, probs))

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "low": self.low,
            "high": self.high,
            "mu": self.mu,
            "sigma": self.sigma,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ReadLengthModel:
        if "kind" not in data and "length" in data:  # shorthand: {"length": 100}
            return cls.fixed(int(data["length"]))
        return cls(
            kind=str(data["kind"]),
            low=int(data["low"]),
            high=int(data["high"]),
            mu=float(data.get("mu", 0.0)),
            sigma=float(data.get("sigma", 0.0)),
        )


@dataclass(frozen=True)
class SimulationConfig:
    """Every parameter of one simulation run."""

    coverage: float
    read_length: ReadLengthModel
    substitution_rate: float = 0.0
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    seed: int = 0

    def __post_init__(self) -> None:
        if self.coverage <= 0:
            raise ValueError("coverage must be positive")
        for name in ("substitution_rate", "insertion_rate", "deletion_rate"):
            rate = getattr(self, name)
            if not 0.0 <= rate <= MAX_ERROR_RATE:
                raise ValueError(f"{name} must lie in [0, {MAX_ERROR_RATE}], got {rate}")

    def n_reads(self, genome_length: int) -> int:
        """Number of reads for a genome: round(coverage * G / mean read length), at least 1."""
        if self.read_length.high > genome_length:
            raise ValueError("reads cannot be longer than the genome")
        return max(1, round(self.coverage * genome_length / self.read_length.mean))

    def nominal_coverage(self, genome_length: int) -> float:
        """c = N * mean_read_length / G for the number of reads actually used."""
        return self.n_reads(genome_length) * self.read_length.mean / genome_length

    def to_dict(self) -> dict[str, Any]:
        return {
            "coverage": self.coverage,
            "read_length": self.read_length.to_dict(),
            "substitution_rate": self.substitution_rate,
            "insertion_rate": self.insertion_rate,
            "deletion_rate": self.deletion_rate,
            "seed": self.seed,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SimulationConfig:
        return cls(
            coverage=float(data["coverage"]),
            read_length=ReadLengthModel.from_dict(data["read_length"]),
            substitution_rate=float(data.get("substitution_rate", 0.0)),
            insertion_rate=float(data.get("insertion_rate", 0.0)),
            deletion_rate=float(data.get("deletion_rate", 0.0)),
            seed=int(data.get("seed", 0)),
        )


@dataclass(frozen=True)
class ReadSet:
    """Simulated reads with their ground truth (structure of arrays, one entry per read)."""

    genome_accession: str
    genome_length: int
    config: SimulationConfig
    sequences: list[str]
    starts: IntArray  # 0-based position of the first reference base of the read window
    spans: IntArray  # number of reference bases the read window covers
    n_substitutions: IntArray
    n_insertions: IntArray
    n_deletions: IntArray

    @property
    def n_reads(self) -> int:
        return len(self.sequences)

    @property
    def ids(self) -> list[str]:
        return [f"read_{i:06d}" for i in range(self.n_reads)]

    @property
    def lengths(self) -> IntArray:
        return np.array([len(s) for s in self.sequences], dtype=np.int64)

    @property
    def realised_coverage(self) -> float:
        """Total reference bases covered by the reads divided by G."""
        return float(self.spans.sum() / self.genome_length)

    def to_fasta(self) -> str:
        """FASTA text; each header carries the ground truth of the read."""
        lines: list[str] = []
        for i, seq in enumerate(self.sequences):
            lines.append(
                f">read_{i:06d} start={self.starts[i]} span={self.spans[i]} "
                f"sub={self.n_substitutions[i]} ins={self.n_insertions[i]} "
                f"del={self.n_deletions[i]}"
            )
            lines.append(seq)
        return "\n".join(lines) + "\n"


def streams(seed: int) -> tuple[np.random.Generator, np.random.Generator]:
    """Independent random streams for read geometry and for errors, both derived from seed."""
    geometry, errors = np.random.SeedSequence(seed).spawn(2)
    return np.random.default_rng(geometry), np.random.default_rng(errors)


def sample_geometry(
    genome_length: int, config: SimulationConfig, rng: np.random.Generator
) -> tuple[IntArray, IntArray]:
    """Draw read lengths and start positions only (no sequences): ``(starts, spans)``."""
    n_reads = config.n_reads(genome_length)
    support, probs = config.read_length.pmf()
    spans = rng.choice(support, size=n_reads, p=probs).astype(np.int64)
    starts = rng.integers(0, genome_length - spans + 1, dtype=np.int64)
    return starts, spans


def _mutate(
    window: NDArray[np.uint8], rng: np.random.Generator, config: SimulationConfig
) -> tuple[NDArray[np.uint8], int, int, int]:
    """Apply the error model to one reference window; returns codes and the three edit counts."""
    n = window.size
    keep = (
        rng.random(n) >= config.deletion_rate
        if config.deletion_rate > 0
        else np.ones(n, dtype=bool)
    )
    bases = window.copy()
    n_sub = 0
    if config.substitution_rate > 0:
        sub = (rng.random(n) < config.substitution_rate) & keep
        n_sub = int(sub.sum())
        if n_sub:
            shift = rng.integers(1, 4, size=n_sub)
            bases[sub] = ((bases[sub].astype(np.int64) + shift) % 4).astype(np.uint8)
    n_ins = 0
    if config.insertion_rate > 0:
        ins = rng.random(n) < config.insertion_rate
        n_ins = int(ins.sum())
        merged = np.empty(2 * n, dtype=np.uint8)
        valid = np.empty(2 * n, dtype=bool)
        merged[0::2], valid[0::2] = bases, keep
        merged[1::2], valid[1::2] = rng.integers(0, 4, size=n).astype(np.uint8), ins
        out = merged[valid]
    else:
        out = bases[keep]
    if out.size == 0:
        raise ValueError("a read lost all its bases to deletions; use a longer read length")
    return out, n_sub, n_ins, int(n - keep.sum())


def simulate_reads(genome: Genome, config: SimulationConfig) -> ReadSet:
    """Simulate a read set from a real genome. Deterministic for a given (genome, config)."""
    geometry_rng, error_rng = streams(config.seed)
    starts, spans = sample_geometry(genome.length, config, geometry_rng)
    codes = _ENCODE[np.frombuffer(genome.sequence.encode("ascii"), dtype=np.uint8)]
    n_reads = starts.size
    has_errors = (config.substitution_rate + config.insertion_rate + config.deletion_rate) > 0
    sequences: list[str] = []
    n_sub = np.zeros(n_reads, dtype=np.int64)
    n_ins = np.zeros(n_reads, dtype=np.int64)
    n_del = np.zeros(n_reads, dtype=np.int64)
    for i in range(n_reads):
        window = codes[starts[i] : starts[i] + spans[i]]
        if has_errors:
            out, n_sub[i], n_ins[i], n_del[i] = _mutate(window, error_rng, config)
        else:
            out = window
        sequences.append(_DECODE[out].tobytes().decode("ascii"))
    return ReadSet(
        genome_accession=genome.accession,
        genome_length=genome.length,
        config=config,
        sequences=sequences,
        starts=starts,
        spans=spans,
        n_substitutions=n_sub,
        n_insertions=n_ins,
        n_deletions=n_del,
    )
