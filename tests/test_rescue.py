"""Tests for the rescue screen design (RESCUE_NEW4_1_48_targets_HQ.m)."""
from pathlib import Path

import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import (
    ANABOLIC_PANEL,
    ANABOLIC_UP_NODES,
    CATABOLIC_DOWN_NODES,
    CATABOLIC_PANEL,
    N_PERTURBATIONS,
    RESCUE_GROUPS,
    enumerate_perturbations,
    run_perturbation,
)

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


@pytest.fixture(scope="module")
def net():
    return load_network(DATA_XLSX)


def test_node_lists_match_manuscript():
    """Methods: six suppression targets (clamped to 0), six activation targets (clamped to 1)."""
    assert CATABOLIC_DOWN_NODES == ("RhoA-E", "PIEZO1", "ROS", "P38", "mTORC1", "NF-κB")
    assert ANABOLIC_UP_NODES == ("SOX9", "PPARγ", "HIF-1α", "NRF2", "SIRT1", "FOXO")


def test_all_rescue_nodes_present_in_network(net):
    targets = CATABOLIC_DOWN_NODES + ANABOLIC_UP_NODES
    assert len(set(targets)) == 12
    for n in targets:
        assert n in net.node_names, f"rescue target {n!r} not in network"
        assert n not in ("Hypo", "NL", "HL")


def test_enumerate_yields_12_singles_and_36_duals():
    perts = list(enumerate_perturbations())
    assert len(perts) == N_PERTURBATIONS == 48
    singles = perts[:12]
    duals = perts[12:]
    assert [p.catabolic_down for p in singles[:6]] == list(CATABOLIC_DOWN_NODES)
    assert [p.anabolic_up for p in singles[6:]] == list(ANABOLIC_UP_NODES)
    assert all(p.anabolic_up and p.catabolic_down for p in duals)
    # MATLAB loop order: activation outer, suppression inner.
    assert duals[0].short_label == "SOX9-RhoA-E"
    assert duals[5].short_label == "SOX9-NF-κB"
    assert duals[-1].short_label == "FOXO-NF-κB"


def test_phenotype_panels_are_network_nodes(net):
    assert len(ANABOLIC_PANEL) == 8 and len(CATABOLIC_PANEL) == 10
    for n in ANABOLIC_PANEL + CATABOLIC_PANEL:
        assert n in net.node_names


def test_rescue_groups_missing_only_non_nodes(net):
    missing = {n for nodes in RESCUE_GROUPS.values() for n in nodes if n not in net.node_names}
    assert missing == {"ATF2", "NutD"}


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
