"""Tests for the falsification module (NP_MT_FALSIFY4_1_39.m)."""
from pathlib import Path

import pytest

from np_mt_rnm.falsification import (
    CONCORDANT,
    CONTRADICTORY,
    FALS_ZERO,
    INDETERMINATE,
    FalsificationRule,
    classify,
    evaluate_rule,
    load_benchmark,
    validate_registry,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

ROOT = Path(__file__).resolve().parents[1]


def test_load_benchmark_parses_39_rules():
    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    assert len(rules) == 39
    up = [r.node for r in rules if r.cls == "anabolic"]
    down = [r.node for r in rules if r.cls == "catabolic"]
    # ANABOLIC_LIST / CATABOLIC_LIST of NP_MT_FALSIFY4_1_39.m.
    assert up == ["ACAN", "COL2A1", "TIMP3", "Bcl2", "SOD2", "SOX9"]
    assert down == [
        "COL1A1", "COL10A1", "MMP3", "MMP13", "ADAMTS4/5",
        "TNF", "IL1β", "IL6", "IL8", "CCL2", "CXCL1", "PGE2", "COX-2",
        "TGFβ", "Wnt3a", "BAX", "CASP3", "CASP9", "HSF1", "MRTF-A",
        "YAP/TAZ", "NF-κB", "β-catenin", "P38", "JNK", "ROS",
        "PIEZO1", "DRP1", "PKN1", "VEGF", "CCN2", "GSK3B", "RUNX2",
    ]
    assert all(r.reference_tag and r.benchmark_ids for r in rules)


def test_registry_maps_onto_network():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    rules = load_benchmark(ROOT / "data" / "falsification_benchmark.csv")
    validate_registry(rules, list(net.node_names))
    with pytest.raises(ValueError):
        validate_registry(rules[1:], list(net.node_names))


@pytest.mark.parametrize("cls, lo, hi, expected", [
    ("anabolic", 0.01, 0.2, CONCORDANT),
    ("anabolic", -0.2, -0.01, CONTRADICTORY),
    ("anabolic", -0.01, 0.2, INDETERMINATE),
    ("catabolic", -0.2, -0.001, CONCORDANT),
    ("catabolic", 0.001, 0.2, CONTRADICTORY),
    ("catabolic", -0.2, 0.0, INDETERMINATE),
])
def test_zero_based_classification(cls, lo, hi, expected):
    """No ±0.02 tolerance: a CI just off zero is decisive."""
    assert FALS_ZERO == 0
    assert classify(cls, lo, hi) == expected


def test_evaluate_rule_nl_higher_passes_for_acan():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    normal = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=20, seed=0, n_jobs=-1)
    hyper = run_replicates(net, REGIME_PRESETS["Hyper"], n_reps=20, seed=0, n_jobs=-1)
    rule = FalsificationRule("ACAN", "anabolic", "Normal>Hyper", "test")
    out = evaluate_rule(rule, normal, hyper, n_boot=500, seed=0)
    assert out.passed and out.outcome == CONCORDANT
    assert out.ci_lower > 0
    assert abs(out.delta - (out.mean_normal - out.mean_hyper)) < 1e-12


def test_evaluate_rule_hl_higher_passes_for_mmp13():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    normal = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=20, seed=0, n_jobs=-1)
    hyper = run_replicates(net, REGIME_PRESETS["Hyper"], n_reps=20, seed=0, n_jobs=-1)
    rule = FalsificationRule("MMP13", "catabolic", "Normal<Hyper", "test")
    out = evaluate_rule(rule, normal, hyper, n_boot=500, seed=0)
    assert out.passed
    assert out.ci_upper < 0
