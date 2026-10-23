"""Operator-selection policies: contracts, behaviour on simulated bandits, forgetting."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genomeasm.bandit import POLICIES, make_policy
from genomeasm.bandit.policies import (
    UCB1,
    EpsilonGreedy,
    Fixed,
    RandomPolicy,
    RouletteWheel,
    RoundRobin,
)


def rng(seed: int = 0) -> np.random.Generator:
    return np.random.default_rng(seed)


def play(policy, means: list[float], pulls: int, generator: np.random.Generator) -> np.ndarray:  # type: ignore[no-untyped-def]
    """Bernoulli bandit: returns how often each arm was chosen."""
    counts = np.zeros(len(means), dtype=int)
    for _ in range(pulls):
        arm = policy.select()
        counts[arm] += 1
        policy.update(arm, float(generator.random() < means[arm]))
    return counts


@pytest.mark.parametrize("name", sorted(POLICIES))
def test_every_policy_obeys_the_contract(name: str) -> None:
    policy = make_policy(name, 4, rng(1))
    for step in range(300):
        arm = policy.select()
        assert 0 <= arm < 4
        p = policy.probabilities()
        assert p.shape == (4,) and abs(p.sum() - 1) < 1e-9 and np.all(p >= 0)
        policy.update(arm, (step % 3) / 2)  # rewards in [0, 1]


def test_factory_rejects_unknown_names_and_bad_parameters() -> None:
    with pytest.raises(ValueError, match="unknown policy"):
        make_policy("magic", 3, rng())
    with pytest.raises(ValueError):
        EpsilonGreedy(3, rng(), epsilon=1.5)
    with pytest.raises(ValueError):
        UCB1(3, rng(), window=0)
    with pytest.raises(ValueError):
        Fixed(3, rng(), arm=3)
    with pytest.raises(ValueError):
        RandomPolicy(0, rng())


def test_random_policy_is_uniform_and_round_robin_cycles() -> None:
    counts = play(RandomPolicy(4, rng(2)), [0.5] * 4, 8000, rng(3))
    assert np.all(np.abs(counts / 8000 - 0.25) < 0.02)
    robin = RoundRobin(3, rng())
    assert [robin.select() for _ in range(7)] == [0, 1, 2, 0, 1, 2, 0]


def test_fixed_policy_always_picks_its_arm() -> None:
    policy = make_policy("fixed", 5, rng(), arm=3)
    assert {policy.select() for _ in range(20)} == {3} and policy.probabilities()[3] == 1.0


@pytest.mark.parametrize("name", ["ucb1", "epsilon_greedy"])
def test_learning_policies_try_every_arm_first(name: str) -> None:
    policy = make_policy(name, 5, rng(4))
    first = []
    for _ in range(5):
        arm = policy.select()
        first.append(arm)
        policy.update(arm, 0.5)
    assert sorted(first) == [0, 1, 2, 3, 4]


@pytest.mark.parametrize("name", ["ucb1", "epsilon_greedy", "roulette"])
def test_learning_policies_favour_the_best_arm(name: str) -> None:
    means = [0.2, 0.5, 0.8]
    counts = play(make_policy(name, 3, rng(5)), means, 4000, rng(6))
    assert int(np.argmax(counts)) == 2 and counts[2] / counts.sum() > 0.5
    assert counts[2] > 2 * counts[0]


def test_random_policy_does_not_favour_anything_but_learners_beat_it() -> None:
    means = [0.1, 0.1, 0.9]
    gain = lambda counts: float(np.dot(counts, means) / counts.sum())  # noqa: E731
    random_gain = gain(play(RandomPolicy(3, rng(7)), means, 3000, rng(8)))
    ucb_gain = gain(play(UCB1(3, rng(7)), means, 3000, rng(8)))
    assert ucb_gain > random_gain + 0.3


def test_sliding_windows_forget_when_the_best_arm_changes() -> None:
    for name in ("ucb1", "epsilon_greedy"):
        policy = make_policy(name, 2, rng(9), window=50)
        generator = rng(10)
        play(policy, [0.9, 0.1], 800, generator)  # arm 0 is best ...
        counts = play(policy, [0.1, 0.9], 400, generator)  # ... then arm 1 becomes best
        assert counts[1] > counts[0], name


def test_roulette_weights_follow_the_rewards_per_segment() -> None:
    policy = RouletteWheel(3, rng(11), segment=10, reaction=0.5)
    for _ in range(10):  # one segment: arm 1 earns 1.0, the others 0
        policy.update(1, 1.0)
        policy.update(0, 0.0)
    assert policy.probabilities()[1] > policy.probabilities()[0]
    assert policy.weights[2] == 1.0  # an unused arm keeps its weight
    for _ in range(200):
        policy.update(0, 0.0)
    assert np.all(policy.weights >= policy.floor)  # no arm disappears


def test_ucb1_index_uses_the_standard_bonus() -> None:
    policy = UCB1(2, rng(), window=100)
    for _ in range(10):
        policy.update(0, 1.0)
    for _ in range(10):
        policy.update(1, 0.0)
    expected = 1.0 + np.sqrt(2 * np.log(20) / 10)  # mean + sqrt(2 ln N / n)
    assert policy._index()[0] == pytest.approx(expected)
    assert policy.select() == 0


@settings(deadline=None, max_examples=40)
@given(seed=st.integers(min_value=0, max_value=10**6), name=st.sampled_from(sorted(POLICIES)))
def test_policies_are_reproducible_per_seed(seed: int, name: str) -> None:
    def run() -> list[int]:
        policy = make_policy(name, 4, rng(seed))
        generator = rng(seed + 1)
        out = []
        for _ in range(100):
            arm = policy.select()
            out.append(arm)
            policy.update(arm, float(generator.random()))
        return out

    assert run() == run()
