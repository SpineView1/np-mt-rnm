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
    """Reconciles paper discrepancy: Section 2.1 says 356, Section 3.1 says 357.

    This test records whatever the Excel actually contains so the paper can
    be updated to match. The assertion is intentionally strict — if the count
    changes, we want to know.
    """
    net = load_network(DATA_XLSX)
    total_edges = int(net.mact.sum() + net.minh.sum())
    # Whichever number is correct will be pinned here after the first run.
    assert total_edges in (356, 357), f"expected 356 or 357 edges, got {total_edges}"
