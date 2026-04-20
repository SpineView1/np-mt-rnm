"""Tests for rescue screen."""
from pathlib import Path

import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import (
    ANABOLIC_UP_NODES,
    CATABOLIC_DOWN_NODES,
    enumerate_perturbations,
    run_perturbation,
)

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


@pytest.fixture(scope="module")
def net():
    return load_network(DATA_XLSX)


def test_node_lists_match_spec():
    """Spec Section 2.5: 5 catabolic targets, 5 anabolic targets.
    Uses the Excel's canonical Greek-letter forms."""
    assert len(CATABOLIC_DOWN_NODES) == 5
    assert len(ANABOLIC_UP_NODES) == 5
    assert set(CATABOLIC_DOWN_NODES) == {"RhoA-E", "PIEZO1", "PI3K-E", "FAK-E", "ROS"}
    assert set(ANABOLIC_UP_NODES) == {"SOX9", "PPARγ", "HIF-1α", "NRF2", "IκBα"}


def test_all_rescue_nodes_present_in_network(net):
    """Sanity: every rescue target must exist in the loaded network."""
    for n in CATABOLIC_DOWN_NODES + ANABOLIC_UP_NODES:
        assert n in net.node_names, f"rescue target {n!r} not in network"


def test_enumerate_yields_35_perturbations():
    perts = list(enumerate_perturbations())
    assert len(perts) == 35


def test_run_perturbation_sox9_shifts_acan_upward(net):
    """SOX9↑ alone should shift ACAN upward vs the unperturbed Hyper baseline."""
    result = run_perturbation(
        net,
        anabolic_up="SOX9",
        catabolic_down=None,
        n_reps=4,
        seed=0,
        n_jobs=1,
    )
    assert result.mean_delta.shape == (len(net.node_names),)
    acan_idx = net.node_names.index("ACAN")
    assert result.mean_delta[acan_idx] > 0
