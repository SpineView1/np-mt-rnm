"""Tests for np_mt_rnm.network — edge-list loader.

Targets: np_mt_rnm/network.py (ports legacy/CreateMatrices_new.m).
"""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


def test_load_network_returns_147_nodes():
    net = load_network(DATA_XLSX)
    assert len(net.node_names) == 147, f"expected 147 nodes, got {len(net.node_names)}"


def test_load_network_mact_minh_shapes():
    net = load_network(DATA_XLSX)
    n = len(net.node_names)
    assert net.mact.shape == (n, n)
    assert net.minh.shape == (n, n)


def test_load_network_edges_are_binary():
    net = load_network(DATA_XLSX)
    assert set(np.unique(net.mact)).issubset({0, 1})
    assert set(np.unique(net.minh)).issubset({0, 1})


def test_load_network_no_self_activation_and_inhibition_overlap():
    """A node cannot be both an activator and inhibitor of the same target."""
    net = load_network(DATA_XLSX)
    overlap = net.mact * net.minh
    assert overlap.sum() == 0, "activator/inhibitor overlap detected"


def test_load_network_stimuli_includes_hypo_nl_hl():
    net = load_network(DATA_XLSX)
    assert "Hypo" in net.stimuli_names
    assert "NL" in net.stimuli_names
    assert "HL" in net.stimuli_names


def test_load_network_edge_count_matches_excel():
    """Pins the edge count of the 2026-09 dataset revision.

    History: the paper's Section 2.1 and abstract say 356 and Section 3.1
    says 357; the pre-revision Excel actually held 357. The 2026-09 revision
    disconnected NutD entirely (4 activation edges removed), leaving 353.
    The paper text needs to be updated to match.
    """
    net = load_network(DATA_XLSX)
    total_edges = int(net.mact.sum() + net.minh.sum())
    assert total_edges == 353, f"expected 353 edges, got {total_edges}"


def test_nutd_is_disconnected():
    """The 2026-09 dataset revision removed every NutD edge.

    Removed: Hypo->NutD, NutD->HIF-1a, NutD->MitD, NutD->ROS (all activation).
    NutD remains in the node list (still 147 nodes) but is now an orphan, so
    it must have no incoming and no outgoing regulation.
    """
    net = load_network(DATA_XLSX)
    i = net.node_names.index("NutD")
    assert net.mact[i, :].sum() == 0 and net.minh[i, :].sum() == 0, "NutD has regulators"
    assert net.mact[:, i].sum() == 0 and net.minh[:, i].sum() == 0, "NutD has targets"


def test_previously_removed_nutd_edges_are_absent():
    """Guards the four specific edges dropped in the 2026-09 revision."""
    net = load_network(DATA_XLSX)
    idx = {n: k for k, n in enumerate(net.node_names)}
    for target, regulator in (
        ("NutD", "Hypo"),
        ("HIF-1\u03b1", "NutD"),
        ("MitD", "NutD"),
        ("ROS", "NutD"),
    ):
        assert net.mact[idx[target], idx[regulator]] == 0, (
            f"edge {regulator} -> {target} should have been removed"
        )
