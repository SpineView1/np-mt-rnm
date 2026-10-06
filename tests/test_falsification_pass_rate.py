"""Falsification benchmark reproduces the paper: 43/45, failing TonEBP and HSF1.

Paper Section 3.4: "Of these, 43 satisfied the predefined criteria ... with
HSF1 and TonEBP as the two exceptions". The Normal and Hyper ensembles are
seeded per replicate exactly as RESCUE_NEW4_1_final.m seeds its baselines
(MATLAB twister, seeds 2002..2101 and 3002..3101).

NP_MT_RNM_FALSIFY4_1.m itself draws rand inside parfor without per-replicate
seeds, so its ensembles (and the exact Δ values in the paper's figure) differ
from run to run even in MATLAB. The pass/fail outcome is robust: the model is
bistable (NRF2/SIRT1/AMPK under Normal, the anabolic attractor under Hyper),
and across 30 independent 100-replicate ensembles 28 gave exactly 43/45 with
TonEBP and HSF1 failing; 2 gave 41/45 (NRF2, SIRT1 also failing).
"""
from pathlib import Path

from np_mt_rnm.falsification import evaluate_benchmark, load_benchmark
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
    assert len(rules) == 45
    outcomes = evaluate_benchmark(rules, normal, hyper, n_boot=10_000, seed=20260420)
    failing = {o.rule.node for o in outcomes if not o.passed}
    assert failing == {"TonEBP", "HSF1"}, failing
