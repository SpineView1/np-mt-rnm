"""Tests for run_single_replicate — one end-to-end ODE solve."""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_single_replicate

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


@pytest.fixture(scope="module")
def net():
    return load_network(DATA_XLSX)


def test_normal_regime_converges_to_steady_state(net):
    result = run_single_replicate(net, regime=REGIME_PRESETS["Normal"], seed=42)
    assert result.converged, f"did not converge; max|dx/dt| = {result.max_abs_derivative:.2e}"
    assert result.max_abs_derivative < 1e-8


def test_clamped_nodes_end_at_clamp_value(net):
    result = run_single_replicate(net, regime=REGIME_PRESETS["Hyper"], seed=123)
    idx_hl = net.node_names.index("HL")
    idx_nl = net.node_names.index("NL")
    idx_hypo = net.node_names.index("Hypo")
    np.testing.assert_allclose(result.x_final[idx_hl], 0.80, atol=1e-10)
    np.testing.assert_allclose(result.x_final[idx_nl], 0.01, atol=1e-10)
    np.testing.assert_allclose(result.x_final[idx_hypo], 0.01, atol=1e-10)


def test_user_clamps_are_enforced(net):
    clamps = {"SOX9": 1.0}
    if "NF-kB" in net.node_names:
        clamps["NF-kB"] = 0.0
    result = run_single_replicate(
        net, regime=REGIME_PRESETS["Hyper"], user_clamps=clamps, seed=7,
    )
    sox9_idx = net.node_names.index("SOX9")
    np.testing.assert_allclose(result.x_final[sox9_idx], 1.0, atol=1e-10)


def test_output_values_in_unit_interval(net):
    result = run_single_replicate(net, regime=REGIME_PRESETS["Normal"], seed=0)
    assert (result.x_final >= -1e-9).all()
    assert (result.x_final <= 1.0 + 1e-9).all()
