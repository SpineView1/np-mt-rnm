"""Tests for falsification module."""
from pathlib import Path

from np_mt_rnm.falsification import (
    FTOL,
    FalsificationRule,
    evaluate_rule,
    load_benchmark,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

ROOT = Path(__file__).resolve().parents[1]


def test_load_benchmark_parses_45_rules():
    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    assert len(rules) == 45, f"expected 45 rules, got {len(rules)}"
    anabolic = [r for r in rules if r.cls == "anabolic"]
    catabolic = [r for r in rules if r.cls == "catabolic"]
    assert len(anabolic) == 17
    assert len(catabolic) == 28


def test_evaluate_rule_anabolic_passes_for_acan():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    normal = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=20, seed=0, n_jobs=-1)
    hyper = run_replicates(net, REGIME_PRESETS["Hyper"], n_reps=20, seed=0, n_jobs=-1)
    rule = FalsificationRule("ACAN", "anabolic", "Normal>Hyper", "test")
    out = evaluate_rule(rule, normal, hyper, n_boot=500, seed=0)
    assert out.passed
    assert out.ci_lower > FTOL


def test_evaluate_rule_catabolic_passes_for_mmp13():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    normal = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=20, seed=0, n_jobs=-1)
    hyper = run_replicates(net, REGIME_PRESETS["Hyper"], n_reps=20, seed=0, n_jobs=-1)
    rule = FalsificationRule("MMP13", "catabolic", "Normal<Hyper", "test")
    out = evaluate_rule(rule, normal, hyper, n_boot=500, seed=0)
    assert out.passed
    assert out.ci_upper < -FTOL
