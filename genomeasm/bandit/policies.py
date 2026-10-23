"""Operator-selection policies for the ALNS (the contribution of this project).

Plain-language story: the ALNS owns several tools for breaking and repairing a solution. Which tool
works best changes during the search, and we do not know in advance. A *policy* chooses the next
tool, sees how well it did (a reward between 0 and 1) and learns. Choosing well while still trying
the other tools from time to time is the *multi-armed bandit* problem: every tool is an arm.

Rewards are in [0, 1]. Because the usefulness of a tool changes over time (the reward distribution
is non-stationary), the learning policies estimate an arm's value from its **last ``window``
rewards** only (a sliding window), not from its whole history.

Policies (``make_policy(name, ...)``):

* ``random``: uniform choice (the baseline the claim is tested against).
* ``roundrobin``: arms in turn.
* ``roulette``: classical ALNS adaptive weights (Ropke and Pisinger 2006; citation TODO-VERIFY). The
  weight of each arm is updated at the end of every ``segment`` iterations to
  ``(1 - reaction) * weight + reaction * (mean reward in the segment)`` (unchanged if the arm was
  not used) and arms are drawn with probability proportional to their weight (minimum weight
  ``floor``, so no arm disappears). An honest strong baseline.
* ``epsilon_greedy``: with probability ``epsilon`` a uniform arm, otherwise the arm with the best
  windowed mean reward (untried arms first).
* ``ucb1``: windowed UCB1 (Auer, Cesa-Bianchi and Fischer 2002; citation TODO-VERIFY): choose the
  maximising ``mean + sqrt(2 ln(N) / n)`` where ``n`` is the number of rewards of the arm in its
  window and ``N`` the total in all windows (untried arms first).
* ``fixed``: always the same arm (a "single operator alone" run).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import deque
from typing import Any

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


class Policy(ABC):
    """Chooses one of ``n_arms`` arms per iteration and learns from the reward in [0, 1]."""

    name = "policy"

    def __init__(self, n_arms: int, rng: np.random.Generator) -> None:
        if n_arms < 1:
            raise ValueError("need at least one arm")
        self.n_arms = n_arms
        self.rng = rng

    @abstractmethod
    def select(self) -> int: ...

    def update(self, arm: int, reward: float) -> None:  # noqa: B027
        """Learn from the reward of the arm that was just used (default: nothing to learn)."""

    @abstractmethod
    def probabilities(self) -> FloatArray:
        """Probability of each arm in the next choice (a one-hot vector for deterministic rules)."""


class RandomPolicy(Policy):
    name = "random"

    def select(self) -> int:
        return int(self.rng.integers(self.n_arms))

    def probabilities(self) -> FloatArray:
        return np.full(self.n_arms, 1.0 / self.n_arms)


class RoundRobin(Policy):
    name = "roundrobin"

    def __init__(self, n_arms: int, rng: np.random.Generator) -> None:
        super().__init__(n_arms, rng)
        self._next = 0

    def select(self) -> int:
        arm = self._next
        self._next = (self._next + 1) % self.n_arms
        return arm

    def probabilities(self) -> FloatArray:
        return np.full(self.n_arms, 1.0 / self.n_arms)


class Fixed(Policy):
    name = "fixed"

    def __init__(self, n_arms: int, rng: np.random.Generator, arm: int = 0) -> None:
        super().__init__(n_arms, rng)
        if not 0 <= arm < n_arms:
            raise ValueError("fixed arm out of range")
        self.arm = arm

    def select(self) -> int:
        return self.arm

    def probabilities(self) -> FloatArray:
        p = np.zeros(self.n_arms)
        p[self.arm] = 1.0
        return p


class RouletteWheel(Policy):
    name = "roulette"

    def __init__(
        self,
        n_arms: int,
        rng: np.random.Generator,
        segment: int = 100,
        reaction: float = 0.5,
        floor: float = 0.01,
    ) -> None:
        super().__init__(n_arms, rng)
        self.segment, self.reaction, self.floor = segment, reaction, floor
        self.weights = np.ones(n_arms)
        self._score = np.zeros(n_arms)
        self._used = np.zeros(n_arms)
        self._iterations = 0

    def select(self) -> int:
        return int(self.rng.choice(self.n_arms, p=self.probabilities()))

    def update(self, arm: int, reward: float) -> None:
        self._score[arm] += reward
        self._used[arm] += 1
        self._iterations += 1
        if self._iterations % self.segment == 0:
            used = self._used > 0
            mean = np.where(used, self._score / np.maximum(self._used, 1), 0.0)
            new = (1 - self.reaction) * self.weights + self.reaction * mean
            self.weights = np.where(used, np.maximum(new, self.floor), self.weights)
            self._score[:] = 0
            self._used[:] = 0

    def probabilities(self) -> FloatArray:
        return self.weights / self.weights.sum()


class _Windowed(Policy):
    """Keeps the last ``window`` rewards of every arm."""

    def __init__(self, n_arms: int, rng: np.random.Generator, window: int = 100) -> None:
        super().__init__(n_arms, rng)
        if window < 1:
            raise ValueError("window must be at least 1")
        self.window = window
        self._rewards = [deque(maxlen=window) for _ in range(n_arms)]  # type: ignore[var-annotated]

    def update(self, arm: int, reward: float) -> None:
        self._rewards[arm].append(reward)

    def counts(self) -> FloatArray:
        return np.array([len(r) for r in self._rewards], dtype=float)

    def means(self) -> FloatArray:
        return np.array([np.mean(r) if r else 0.0 for r in self._rewards])


class EpsilonGreedy(_Windowed):
    name = "epsilon_greedy"

    def __init__(
        self, n_arms: int, rng: np.random.Generator, window: int = 100, epsilon: float = 0.1
    ) -> None:
        super().__init__(n_arms, rng, window)
        if not 0.0 <= epsilon <= 1.0:
            raise ValueError("epsilon must lie in [0, 1]")
        self.epsilon = epsilon

    def select(self) -> int:
        counts = self.counts()
        if (counts == 0).any():
            return int(np.flatnonzero(counts == 0)[0])
        if self.rng.random() < self.epsilon:
            return int(self.rng.integers(self.n_arms))
        means = self.means()
        return int(self.rng.choice(np.flatnonzero(means == np.max(means))))

    def probabilities(self) -> FloatArray:
        counts = self.counts()
        if (counts == 0).any():
            p = (counts == 0).astype(float)
            return p / p.sum()
        means = self.means()
        best = means == np.max(means)
        return self.epsilon / self.n_arms + (1 - self.epsilon) * best / best.sum()


class UCB1(_Windowed):
    name = "ucb1"

    def __init__(
        self, n_arms: int, rng: np.random.Generator, window: int = 100, c: float = 1.0
    ) -> None:
        super().__init__(n_arms, rng, window)
        self.c = c

    def _index(self) -> FloatArray:
        counts = self.counts()
        total = counts.sum()
        return self.means() + self.c * np.sqrt(
            2.0 * np.log(max(total, 2.0)) / np.maximum(counts, 1)
        )

    def select(self) -> int:
        counts = self.counts()
        if (counts == 0).any():
            return int(np.flatnonzero(counts == 0)[0])
        return int(np.argmax(self._index()))

    def probabilities(self) -> FloatArray:
        counts = self.counts()
        p = np.zeros(self.n_arms)
        p[
            int(np.flatnonzero(counts == 0)[0])
            if (counts == 0).any()
            else int(np.argmax(self._index()))
        ] = 1.0
        return p


POLICIES: dict[str, type[Policy]] = {
    "random": RandomPolicy,
    "roundrobin": RoundRobin,
    "roulette": RouletteWheel,
    "epsilon_greedy": EpsilonGreedy,
    "ucb1": UCB1,
    "fixed": Fixed,
}


def make_policy(name: str, n_arms: int, rng: np.random.Generator, **params: Any) -> Policy:
    """Create a policy by name; ``params`` are its keyword parameters (window, epsilon, ...)."""
    try:
        return POLICIES[name](n_arms, rng, **params)
    except KeyError:
        raise ValueError(f"unknown policy {name!r}; choose from {sorted(POLICIES)}") from None
