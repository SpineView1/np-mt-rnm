"""End-to-end invariant: running the full benchmark at 100 reps yields >=95% pass rate.

Paper claims 95.6% (43/45). This test catches regressions where the port
drops below that threshold — a 5k bootstrap is used here to keep CI
reasonable (~2–4 min), sufficient to catch substantial changes.
"""
from pathlib import Path

from np_mt_rnm.falsification import (
    evaluate_benchmark,
    load_benchmark,
    pass_rate,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

ROOT = Path(__file__).resolve().parents[1]


def test_falsification_pass_rate_at_least_95_percent():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    normal = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=100, seed=20260420, n_jobs=-1)
    hyper = run_replicates(
        net, REGIME_PRESETS["Hyper"], n_reps=100, seed=20260420 + 10000, n_jobs=-1
    )
    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    outcomes = evaluate_benchmark(rules, normal, hyper, n_boot=5000, seed=20260420)
    pr = pass_rate(outcomes)
    assert pr >= 0.95, f"pass rate {pr:.3f} below paper's 95% claim"
