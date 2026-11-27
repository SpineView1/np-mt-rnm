"""Falsification against 39 literature-derived benchmarks (manuscript Figure 5).

Ports the falsification part of legacy/NP_MT_FALSIFY4_1_39.m: Δ = mean(NL) −
mean(HL), 95% percentile bootstrap CI (10,000 resamples), zero-based
directional criteria. Writes results/tables/falsification_rules.csv; the
figure panels are drawn by scripts/run_paper_figures.py.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from np_mt_rnm.falsification import (
    CLASS_LABELS,
    FALS_NBOOT,
    evaluate_benchmark,
    load_benchmark,
    outcome_counts,
    pass_rate,
    validate_registry,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates

ROOT = Path(__file__).resolve().parents[1]
N_REPS = 100
SEED = 20260420


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    validate_registry(rules, list(net.node_names))

    print("[falsify] running Normal ensemble ...")
    normal = run_replicates(
        net, REGIME_PRESETS["Normal"], n_reps=N_REPS, seed=REGIME_SEEDS["Normal"], n_jobs=-1
    )
    print("[falsify] running Hyper ensemble ...")
    hyper = run_replicates(
        net, REGIME_PRESETS["Hyper"], n_reps=N_REPS, seed=REGIME_SEEDS["Hyper"], n_jobs=-1
    )

    print(f"[falsify] evaluating {len(rules)} benchmarks with {FALS_NBOOT} bootstrap resamples ...")
    outcomes = evaluate_benchmark(rules, normal, hyper, n_boot=FALS_NBOOT, seed=SEED)

    counts = outcome_counts(outcomes)
    passed = sum(o.passed for o in outcomes)
    print(f"[falsify] concordant {passed}/{len(outcomes)} ({100 * pass_rate(outcomes):.1f}%); "
          f"contradictory {counts['Contradictory']}; indeterminate {counts['Indeterminate']}")
    for o in outcomes:
        if not o.passed:
            print(f"  {o.rule.node} ({CLASS_LABELS[o.rule.cls]}): {o.outcome}, Δ={o.delta:+.4f}, "
                  f"CI=[{o.ci_lower:+.4f}, {o.ci_upper:+.4f}]")

    df = pd.DataFrame([
        {
            "node": o.rule.node,
            "class": o.rule.cls,
            "direction": CLASS_LABELS[o.rule.cls],
            "literature_direction": o.rule.literature_direction,
            "expected": o.rule.expected,
            "benchmark_ids": o.rule.benchmark_ids,
            "reference": o.rule.reference_tag,
            "mean_normal": o.mean_normal,
            "mean_hyper": o.mean_hyper,
            "delta": o.delta,
            "ci_lower": o.ci_lower,
            "ci_upper": o.ci_upper,
            "outcome": o.outcome,
            "passed": o.passed,
        }
        for o in outcomes
    ])
    df.to_csv(ROOT / "results" / "tables" / "falsification_rules.csv", index=False)


if __name__ == "__main__":
    main()
