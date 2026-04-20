"""Tests for simulation.py helpers: regime presets, clamp construction."""
import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import (
    REGIME_PRESETS,
    build_clamps,
)
from pathlib import Path

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


def test_regime_presets_match_spec():
    """Per spec Section 6: Hypo = (0.20, 0.01, 0.01); Normal = (0.01, 0.80, 0.01);
    Hyper = (0.01, 0.01, 0.80)."""
    assert REGIME_PRESETS["Hypo"] == {"Hypo": 0.20, "NL": 0.01, "HL": 0.01}
    assert REGIME_PRESETS["Normal"] == {"Hypo": 0.01, "NL": 0.80, "HL": 0.01}
    assert REGIME_PRESETS["Hyper"] == {"Hypo": 0.01, "NL": 0.01, "HL": 0.80}


def test_build_clamps_marks_regime_inputs():
    net = load_network(DATA_XLSX)
    clamped_mask, x_clamp = build_clamps(
        net, regime=REGIME_PRESETS["Normal"], user_clamps=None
    )
    for name in ("Hypo", "NL", "HL"):
        idx = net.node_names.index(name)
        assert clamped_mask[idx], f"{name} should be clamped"
    nl_idx = net.node_names.index("NL")
    assert x_clamp[nl_idx] == 0.80


def test_build_clamps_honors_user_clamps():
    net = load_network(DATA_XLSX)
    clamped_mask, x_clamp = build_clamps(
        net, regime=REGIME_PRESETS["Hyper"], user_clamps={"SOX9": 1.0}
    )
    sox9_idx = net.node_names.index("SOX9")
    assert clamped_mask[sox9_idx]
    assert x_clamp[sox9_idx] == 1.0


def test_build_clamps_unknown_user_node_raises():
    net = load_network(DATA_XLSX)
    with pytest.raises(KeyError):
        build_clamps(net, regime=REGIME_PRESETS["Normal"], user_clamps={"NOT_A_NODE": 1})
