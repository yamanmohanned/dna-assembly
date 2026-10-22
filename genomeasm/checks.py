"""Small helpers shared by the validation suites: recording checks and basic statistics."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

Checks = list[dict[str, Any]]


def plain(value: Any) -> Any:
    """Turn numpy scalars into plain Python numbers so the result can be written as JSON."""
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def add_check(
    checks: Checks, check_id: str, description: str, passed: bool, **numbers: Any
) -> None:
    """Append one check (id, description, verdict and the numbers behind it) to the list."""
    checks.append(
        {
            "id": check_id,
            "description": description,
            "passed": bool(passed),
            **{k: plain(v) for k, v in numbers.items()},
        }
    )


def mean_and_se(values: Sequence[float] | np.ndarray) -> tuple[float, float]:
    """Mean and standard error (sample standard deviation divided by the square root of n)."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.size < 2:
        return float(arr.mean()), 0.0
    return float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(arr.size))
