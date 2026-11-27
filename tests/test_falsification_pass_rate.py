"""Falsification reproduces the manuscript: 37/39 concordant, HSF1 and CCN2 contradictory.

Manuscript Results: "Of the 39 experimental benchmarks, 37 (94.9%) showed
RNM-predicted effects concordant with their literature-derived directionality
... The remaining two benchmarks, HSF1 and CCN2, were contradictory (5.1%), and
none was indeterminate." The Normal and Hyper ensembles are seeded per
replicate exactly as RESCUE_NEW4_1_48_targets_HQ.m seeds its baselines (MATLAB
twister, seeds 2002..2101 and 3002..3101).

NP_MT_FALSIFY4_1_39.m itself draws rand inside parfor without per-replicate
seeds, so its ensembles (and the exact Δ values in the manuscript's
Supplementary Table) differ from run to run even in MATLAB; the outcome does
not.
"""
from pathlib import Path

from np_mt_rnm.falsification import (
    evaluate_benchmark,
    load_benchmark,
    outcome_counts,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates

ROOT = Path(__file__).resolve().parents[1]


def test_falsification_matches_paper():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    normal = run_replicates(
        net, REGIME_PRESETS["Normal"], n_reps=100, seed=REGIME_SEEDS["Normal"], n_jobs=-1
    )
    hyper = run_replicates(
        net, REGIME_PRESETS["Hyper"], n_reps=100, seed=REGIME_SEEDS["Hyper"], n_jobs=-1
    )
    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    assert len(rules) == 39
    outcomes = evaluate_benchmark(rules, normal, hyper, n_boot=10_000, seed=20260420)
    assert outcome_counts(outcomes) == {"Concordant": 37, "Contradictory": 2, "Indeterminate": 0}
    failing = {o.rule.node for o in outcomes if not o.passed}
    assert failing == {"HSF1", "CCN2"}, failing
    up = [o for o in outcomes if o.rule.cls == "anabolic"]
    assert all(o.passed for o in up)  # 6/6
