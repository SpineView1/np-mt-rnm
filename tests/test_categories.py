"""Tests for categories.py.

Ensures every node in the edge list has at least one category, and that
categories are consistent with the MATLAB group_categories source.
"""
from pathlib import Path

import pytest

from np_mt_rnm.categories import (
    CATEGORY_LABELS,
    CATEGORY_ORDER,
    NODE_CATEGORIES,
    categories_for,
    get_primary_category,
    nodes_in_category,
)
from np_mt_rnm.network import load_network

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


# NutD was disconnected from the network in the 2026-09 dataset revision and
# dropped from group_categories in NP_MT_RNM_FALSIFY4_1.m at the same time, so
# it is deliberately in no category. (RESCUE_NEW4_1_final.m still lists it under
# "Mechanical stimuli & its receptors" — flagged to Zerihun as an inconsistency.)
UNCATEGORISED_BY_DESIGN = {"NutD"}


def test_every_node_has_at_least_one_category():
    """Every connected node in MT_PRIMARY4_1.xlsx must have a category.

    Missing nodes indicate either the Excel has nodes not listed in the MATLAB
    group_categories, OR the transcription from MATLAB to Python dropped some.
    Either way — surface it, don't hide it. The only permitted exceptions are
    nodes the biologist explicitly removed from group_categories.
    """
    net = load_network(DATA_XLSX)
    missing = [
        n for n in net.node_names
        if n not in NODE_CATEGORIES and n not in UNCATEGORISED_BY_DESIGN
    ]
    assert not missing, (
        f"{len(missing)} of 147 nodes have no category assignment:\n"
        + "\n".join(f"  - {n}" for n in missing)
    )


def test_uncategorised_nodes_are_disconnected():
    """A node may only be uncategorised if it has no edges at all."""
    net = load_network(DATA_XLSX)
    for name in UNCATEGORISED_BY_DESIGN:
        i = net.node_names.index(name)
        degree = int(
            net.mact[i, :].sum() + net.minh[i, :].sum()
            + net.mact[:, i].sum() + net.minh[:, i].sum()
        )
        assert degree == 0, f"{name} is uncategorised but still has {degree} edges"


def test_category_order_covers_every_declared_category():
    assert set(CATEGORY_ORDER) == set(CATEGORY_LABELS.keys())


def test_nodes_in_category_returns_lists_from_matlab_source():
    """Sanity spot-checks that the transcription matches MATLAB."""
    # Spot-check a few lists that are short and unique:
    ecm = nodes_in_category("ecm_matrix")
    assert ecm == ["COL2A1", "COL1A1", "COL10A1", "ACAN", "TIMP3"]

    mech = nodes_in_category("mechanosensor")
    assert "PIEZO1" in mech
    assert "TRPV4" in mech
    assert "Hypo" in mech and "NL" in mech and "HL" in mech


def test_primary_category_is_first_per_category_order():
    """HO-1 appears in both 'metabolic' and 'oxidative_proteostasis'.
    CATEGORY_ORDER puts metabolic before oxidative_proteostasis, so
    primary for HO-1 should be 'metabolic'."""
    assert get_primary_category("HO-1") == "metabolic"
    assert "oxidative_proteostasis" in categories_for("HO-1")


def test_unknown_node_raises():
    with pytest.raises(KeyError):
        get_primary_category("NOT_A_NODE")


def test_unknown_category_raises():
    with pytest.raises(KeyError):
        nodes_in_category("not_a_category")


def test_other_category_contains_only_matlab_unassigned_nodes():
    """The 'other' category is the catch-all for nodes not in any MATLAB
    group_categories. This preserves the MATLAB source's silence — we don't
    invent biological assignments."""
    other = nodes_in_category("other")
    expected = {
        "ANXA1", "BAD", "DAG", "IKK", "iNOS", "IκBα",
        "LIMK1", "MAP3K", "SMIT1", "SOX5", "SOX6", "TAUT",
        "β-catenin",
    }
    assert set(other) == expected
