"""Reproduce paper Figure 9 + falsification pass-rate."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from np_mt_rnm.falsification import (
    evaluate_benchmark,
    load_benchmark,
    pass_rate,
)
from np_mt_rnm.figures import plot_fig9_forest
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

ROOT = Path(__file__).resolve().parents[1]
N_REPS = 100
N_BOOT = 10_000
SEED = 20260420


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    print("[falsify] running Normal ensemble ...")
    normal = run_replicates(
        net, REGIME_PRESETS["Normal"], n_reps=N_REPS, seed=SEED, n_jobs=-1
    )
    print("[falsify] running Hyper ensemble ...")
    hyper = run_replicates(
        net,
        REGIME_PRESETS["Hyper"],
        n_reps=N_REPS,
        seed=SEED + 10000,
        n_jobs=-1,
    )

    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    print(f"[falsify] evaluating {len(rules)} rules with {N_BOOT} bootstrap resamples ...")
    outcomes = evaluate_benchmark(rules, normal, hyper, n_boot=N_BOOT, seed=SEED)

    pr = pass_rate(outcomes)
    passed = sum(o.passed for o in outcomes)
    print(f"[falsify] pass rate: {pr:.3f} ({passed}/{len(outcomes)})")
    if pr < 0.95:
        failing = [o for o in outcomes if not o.passed]
        print("[falsify] NOT PASSED (< 95%). Failing rules:")
        for o in failing:
            print(f"  {o.rule.node} ({o.rule.cls}): Δ={o.delta:.3f}, CI=[{o.ci_lower:.3f}, {o.ci_upper:.3f}]")

    plot_fig9_forest(outcomes, ROOT / "results" / "figures" / "fig9_falsification_forest.png")

    df = pd.DataFrame([
        {
            "node": o.rule.node,
            "class": o.rule.cls,
            "expected": o.rule.expected,
            "reference": o.rule.reference_tag,
            "delta": o.delta,
            "ci_lower": o.ci_lower,
            "ci_upper": o.ci_upper,
            "passed": o.passed,
        }
        for o in outcomes
    ])
    df.to_csv(ROOT / "results" / "tables" / "falsification_rules.csv", index=False)


if __name__ == "__main__":
    main()
