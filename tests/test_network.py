"""Tests for np_mt_rnm.network — edge-list loader.

Targets: np_mt_rnm/network.py (ports legacy/CreateMatrices_new.m).
"""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


def test_load_network_returns_146_nodes():
    net = load_network(DATA_XLSX)
    assert len(net.node_names) == 146, f"expected 146 nodes, got {len(net.node_names)}"


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
    """The network is the 353-interaction corpus (Supplementary Data 1).

    353 directed interactions (277 activation + 76 inhibition) between 146
    nodes. An earlier 357-edge version also had four NutD interactions
    (Hypo->NutD, NutD->HIF-1a, NutD->MitD, NutD->ROS); they are not part of the
    curated corpus, so they and the NutD node were removed.
    """
    net = load_network(DATA_XLSX)
    assert int(net.mact.sum()) == 277
    assert int(net.minh.sum()) == 76
    assert int(net.mact.sum() + net.minh.sum()) == 353


def test_nutd_removed():
    """NutD had no interactions left in the corpus, so it is not a node."""
    net = load_network(DATA_XLSX)
    assert "NutD" not in net.node_names
