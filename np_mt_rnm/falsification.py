"""Falsification benchmark. Ports legacy/NP_MT_FALSIFY4_1_39.m (falsification section).

Per the manuscript (Methods, "Model falsification using independent NP
mechanobiology experiments"):

  - 39 literature-derived benchmarks: 6 markers with higher expression/activity
    under NL-like loading (literature direction "+") and 33 with higher
    expression/activity under HL-like loading ("-"). See
    data/falsification_benchmark.csv.
  - For each benchmark node, Δ = mean(NL) − mean(HL) across 100 replicates.
  - 95% percentile bootstrap CI of Δ, N_boot = 10,000, resampling the NL and HL
    replicate outcomes independently.
  - Zero-based directional criteria (no tolerance band):
      "+" benchmark: concordant if CI_lower > 0, contradictory if CI_upper < 0;
      "-" benchmark: concordant if CI_upper < 0, contradictory if CI_lower > 0;
      otherwise (CI includes zero) indeterminate.

The legacy internal labels "anabolic" / "catabolic" are kept as the ``cls``
values of the registry; they mean an upward ("+", NL-higher) or downward
("-", HL-higher) benchmark response, not that every marker is biologically
anabolic or catabolic.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.simulation import ReplicateEnsemble
from np_mt_rnm.statistics import bootstrap_diff_ci

# FALS_ZERO in NP_MT_FALSIFY4_1_39.m: zero is the only reference.
FALS_ZERO = 0.0
FALS_ALPHA = 0.05
FALS_NBOOT = 10_000

N_NL_HIGHER = 6
N_HL_HIGHER = 33
N_BENCHMARKS = N_NL_HIGHER + N_HL_HIGHER  # 39

CONCORDANT = "Concordant"
CONTRADICTORY = "Contradictory"
INDETERMINATE = "Indeterminate"

# Display labels for the internal classes.
CLASS_LABELS = {"anabolic": "NL-higher", "catabolic": "HL-higher"}


@dataclass(frozen=True)
class FalsificationRule:
    node: str
    cls: str                # "anabolic" (literature "+") or "catabolic" ("-")
    expected: str           # "Normal>Hyper" or "Normal<Hyper"
    reference_tag: str      # ';'-separated BibTeX keys of the supporting studies
    benchmark_ids: str = ""  # ';'-separated benchmark entry IDs (B01..B51)

    @property
    def literature_direction(self) -> str:
        return "+" if self.cls == "anabolic" else "-"


@dataclass(frozen=True)
class FalsificationOutcome:
    rule: FalsificationRule
    delta: float
    ci_lower: float
    ci_upper: float
    mean_normal: float
    mean_hyper: float
    outcome: str            # Concordant / Contradictory / Indeterminate

    @property
    def passed(self) -> bool:
        return self.outcome == CONCORDANT


def load_benchmark(csv_path: Path | str) -> list[FalsificationRule]:
    df = pd.read_csv(csv_path, keep_default_na=False)
    required = {"node", "class", "expected_direction", "reference_tag"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"benchmark CSV missing columns: {missing}")
    rules = [
        FalsificationRule(
            node=str(row["node"]),
            cls=str(row["class"]),
            expected=str(row["expected_direction"]),
            reference_tag=str(row["reference_tag"]),
            benchmark_ids=str(row.get("benchmark_ids", "")),
        )
        for _, row in df.iterrows()
    ]
    nodes = [r.node for r in rules]
    if len(set(nodes)) != len(nodes):
        raise ValueError("benchmark registry lists a node more than once")
    for r in rules:
        if r.cls not in CLASS_LABELS:
            raise ValueError(f"rule class must be 'anabolic' or 'catabolic', got {r.cls!r}")
    return rules


def validate_registry(rules: list[FalsificationRule], node_names: list[str]) -> None:
    """validate_benchmark_registry in NP_MT_FALSIFY4_1_39.m: fail before simulating.

    Requires exactly 6 NL-higher and 33 HL-higher distinct benchmarks, each an
    existing network node.
    """
    n_up = sum(r.cls == "anabolic" for r in rules)
    n_down = sum(r.cls == "catabolic" for r in rules)
    if (n_up, n_down) != (N_NL_HIGHER, N_HL_HIGHER):
        raise ValueError(
            f"registry must contain {N_NL_HIGHER} NL-higher and {N_HL_HIGHER} "
            f"HL-higher markers, found {n_up} and {n_down}"
        )
    absent = [r.node for r in rules if r.node not in node_names]
    if absent:
        raise ValueError(f"benchmark nodes not in the network: {absent}")


def classify(cls: str, ci_lower: float, ci_upper: float) -> str:
    """Zero-based CI classification of one benchmark."""
    if ci_lower > FALS_ZERO:
        observed = "+"
    elif ci_upper < FALS_ZERO:
        observed = "-"
    else:
        return INDETERMINATE
    expected = "+" if cls == "anabolic" else "-"
    return CONCORDANT if observed == expected else CONTRADICTORY


def evaluate_rule(
    rule: FalsificationRule,
    normal_ens: ReplicateEnsemble,
    hyper_ens: ReplicateEnsemble,
    n_boot: int = FALS_NBOOT,
    seed: int = 0,
) -> FalsificationOutcome:
    if rule.node not in normal_ens.node_names:
        raise KeyError(f"rule node {rule.node!r} missing from ensemble")
    if rule.cls not in CLASS_LABELS:
        raise ValueError(f"rule class must be 'anabolic' or 'catabolic', got {rule.cls!r}")
    j = normal_ens.node_names.index(rule.node)
    x = normal_ens.steady_states[:, j]
    y = hyper_ens.steady_states[:, j]
    rng = np.random.default_rng(seed)
    delta, lo, hi = bootstrap_diff_ci(x, y, n_boot=n_boot, alpha=FALS_ALPHA, rng=rng)
    return FalsificationOutcome(
        rule=rule,
        delta=delta,
        ci_lower=lo,
        ci_upper=hi,
        mean_normal=float(x.mean()),
        mean_hyper=float(y.mean()),
        outcome=classify(rule.cls, lo, hi),
    )


def evaluate_benchmark(
    rules: list[FalsificationRule],
    normal_ens: ReplicateEnsemble,
    hyper_ens: ReplicateEnsemble,
    n_boot: int = FALS_NBOOT,
    seed: int = 0,
) -> list[FalsificationOutcome]:
    return [
        evaluate_rule(r, normal_ens, hyper_ens, n_boot=n_boot, seed=seed + i)
        for i, r in enumerate(rules)
    ]


def pass_rate(outcomes: list[FalsificationOutcome]) -> float:
    return sum(o.passed for o in outcomes) / len(outcomes)


def outcome_counts(outcomes: list[FalsificationOutcome]) -> dict[str, int]:
    """Supported / contradicted / inconclusive counts (Figure B key)."""
    return {
        CONCORDANT: sum(o.outcome == CONCORDANT for o in outcomes),
        CONTRADICTORY: sum(o.outcome == CONTRADICTORY for o in outcomes),
        INDETERMINATE: sum(o.outcome == INDETERMINATE for o in outcomes),
    }
