"""Tests for transitions.py — constrained Hypo→Normal and Normal→Hyper paths."""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.transitions import (
    HYPO_TO_NORMAL_PATH,
    NORMAL_TO_HYPER_PATH,
    run_transition_path,
)

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


@pytest.fixture(scope="module")
def net():
    return load_network(DATA_XLSX)


def test_hypo_to_normal_path_is_six_steps():
    assert len(HYPO_TO_NORMAL_PATH) == 6
    assert HYPO_TO_NORMAL_PATH[0]["Hypo"] > HYPO_TO_NORMAL_PATH[0]["NL"]
    assert HYPO_TO_NORMAL_PATH[-1]["NL"] > HYPO_TO_NORMAL_PATH[-1]["Hypo"]


def test_normal_to_hyper_path_is_six_steps():
    assert len(NORMAL_TO_HYPER_PATH) == 6
    assert NORMAL_TO_HYPER_PATH[0]["NL"] > NORMAL_TO_HYPER_PATH[0]["HL"]
    assert NORMAL_TO_HYPER_PATH[-1]["HL"] > NORMAL_TO_HYPER_PATH[-1]["NL"]


def test_run_transition_path_shape(net):
    result = run_transition_path(
        net, HYPO_TO_NORMAL_PATH, n_reps=4, seed=0, n_jobs=1
    )
    assert result.means.shape == (6, len(net.node_names))
    assert result.stds.shape == (6, len(net.node_names))


def test_hypo_to_normal_anabolic_endpoint_greater_than_start(net):
    """Per Fig 6/7, ACAN/COL2A1 should finish higher than they start along Hypo→Normal."""
    result = run_transition_path(net, HYPO_TO_NORMAL_PATH, n_reps=8, seed=0, n_jobs=-1)
    acan_idx = net.node_names.index("ACAN")
    col2a1_idx = net.node_names.index("COL2A1")
    assert result.means[-1, acan_idx] > result.means[0, acan_idx]
    assert result.means[-1, col2a1_idx] > result.means[0, col2a1_idx]


def test_run_trajectory_shape_and_switch():
    """Trajectory has expected shape and reflects regime switch at t_switch."""
    from np_mt_rnm.simulation import REGIME_PRESETS
    from np_mt_rnm.transitions import run_trajectory

    net = load_network(DATA_XLSX)
    result = run_trajectory(
        net,
        regime_a=REGIME_PRESETS["Normal"],
        regime_b=REGIME_PRESETS["Hyper"],
        t_switch=10.0,
        t_end=30.0,
        n_t=31,
        n_reps=2,
        seed=42,
        n_jobs=1,
    )
    assert result.t.shape == (31,)
    assert result.mean.shape == (31, len(net.node_names))
    assert result.std.shape == (31, len(net.node_names))
    assert result.t[0] == 0.0
    assert result.t[-1] == 30.0
    # t_switch=10 in a 31-pt grid 0..30 (step 1.0) lands at index 10.
    assert np.isclose(result.t[10], 10.0)
