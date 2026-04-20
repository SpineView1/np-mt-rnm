"""Falsification benchmark. Ports legacy/NP_MT_RNM_FALSIFY4_1.m (falsification section).

Per paper Section 2.4:
  - For each benchmark node, compute Δ = mean(Normal) − mean(Hyper) across replicates.
  - Bootstrap 95% CI for Δ with 10k resamples (FALS_NBOOT).
  - Anabolic rule passes if CI_lower > +FTOL (FTOL = 0.02).
  - Catabolic rule passes if CI_upper < -FTOL.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.simulation import ReplicateEnsemble
from np_mt_rnm.statistics import bootstrap_diff_ci

FTOL = 0.02


@dataclass(frozen=True)
class FalsificationRule:
    node: str
    cls: str                # "anabolic" or "catabolic"
    expected: str           # "Normal>Hyper" or "Normal<Hyper"
    reference_tag: str


@dataclass(frozen=True)
class FalsificationOutcome:
    rule: FalsificationRule
    delta: float
    ci_lower: float
    ci_upper: float
    passed: bool


def load_benchmark(csv_path: Path | str) -> list[FalsificationRule]:
    df = pd.read_csv(csv_path)
    required = {"node", "class", "expected_direction", "reference_tag"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"benchmark CSV missing columns: {missing}")
    return [
        FalsificationRule(
            node=str(row["node"]),
            cls=str(row["class"]),
            expected=str(row["expected_direction"]),
            reference_tag=str(row["reference_tag"]),
        )
        for _, row in df.iterrows()
    ]


def evaluate_rule(
    rule: FalsificationRule,
    normal_ens: ReplicateEnsemble,
    hyper_ens: ReplicateEnsemble,
    n_boot: int = 10_000,
    seed: int = 0,
) -> FalsificationOutcome:
    if rule.node not in normal_ens.node_names:
        raise KeyError(f"rule node {rule.node!r} missing from ensemble")
    j = normal_ens.node_names.index(rule.node)
    x = normal_ens.steady_states[:, j]
    y = hyper_ens.steady_states[:, j]
    rng = np.random.default_rng(seed)
    delta, lo, hi = bootstrap_diff_ci(x, y, n_boot=n_boot, alpha=0.05, rng=rng)
    if rule.cls == "anabolic":
        passed = lo > FTOL
    elif rule.cls == "catabolic":
        passed = hi < -FTOL
    else:
        raise ValueError(f"rule class must be 'anabolic' or 'catabolic', got {rule.cls!r}")
    return FalsificationOutcome(rule=rule, delta=delta, ci_lower=lo, ci_upper=hi, passed=passed)


def evaluate_benchmark(
    rules: list[FalsificationRule],
    normal_ens: ReplicateEnsemble,
    hyper_ens: ReplicateEnsemble,
    n_boot: int = 10_000,
    seed: int = 0,
) -> list[FalsificationOutcome]:
    return [
        evaluate_rule(r, normal_ens, hyper_ens, n_boot=n_boot, seed=seed + i)
        for i, r in enumerate(rules)
    ]


def pass_rate(outcomes: list[FalsificationOutcome]) -> float:
    return sum(o.passed for o in outcomes) / len(outcomes)
