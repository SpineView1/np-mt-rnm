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
    """Pins the network the manuscript's figures were generated with.

    357 directed edges (281 activation + 76 inhibition). Every rescue
    percentage printed in the manuscript (text and Top-20 figures) is
    reproduced to the printed decimal with this network and not with the
    353-edge 2026-09 revision, which removed NutD's four edges; see
    tests/test_paper_reproduction.py.
    """
    net = load_network(DATA_XLSX)
    assert int(net.mact.sum()) == 281
    assert int(net.minh.sum()) == 76
    assert int(net.mact.sum() + net.minh.sum()) == 357


def test_nutd_edges_present():
    """The four NutD edges dropped in the 353-edge revision are present."""
    net = load_network(DATA_XLSX)
    idx = {n: k for k, n in enumerate(net.node_names)}
    for target, regulator in (
        ("NutD", "Hypo"),
        ("HIF-1\u03b1", "NutD"),
        ("MitD", "NutD"),
        ("ROS", "NutD"),
    ):
        assert net.mact[idx[target], idx[regulator]] == 1, (
            f"edge {regulator} -> {target} missing"
        )
