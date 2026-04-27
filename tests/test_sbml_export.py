"""SBML export sanity tests."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.sbml_export import _build_omega_formula, _sanitize_id, export_sbml

DATA = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


def test_sanitize_id_handles_greek_and_specials():
    assert _sanitize_id("HIF-1α") == "HIF_1alpha"
    assert _sanitize_id("NF-κB") == "NF_kappaB"
    assert _sanitize_id("β-catenin") == "beta_catenin"
    assert _sanitize_id("PPARγ") == "PPARgamma"
    assert _sanitize_id("ADAMTS4/5") == "ADAMTS4_5"
    # Leading digit gets prefixed underscore
    assert _sanitize_id("3PG").startswith("_")


def test_omega_formula_three_cases():
    only_act = _build_omega_formula(["A", "B"], [])
    assert "A + B" in only_act
    only_inh = _build_omega_formula([], ["X"])
    assert only_inh.startswith("(1 - ")
    both = _build_omega_formula(["A"], ["X"])
    assert "1 -" in both
    none = _build_omega_formula([], [])
    assert none == "0"


def test_export_produces_valid_sbml(tmp_path):
    libsbml = pytest.importorskip("libsbml")
    out = tmp_path / "model.xml"
    net = load_network(DATA)
    ok = export_sbml(net, out)
    assert ok
    assert out.exists() and out.stat().st_size > 1000

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(str(out))
    assert doc.getNumErrors() == 0
    m = doc.getModel()
    assert m.getNumSpecies() == 147
    # Boundary species are exactly the 3 mechanical loading inputs.
    boundary = sorted(
        m.getSpecies(i).getId()
        for i in range(m.getNumSpecies())
        if m.getSpecies(i).getBoundaryCondition()
    )
    assert boundary == sorted(["Hypo", "NL", "HL"])
    # Rate rules: one per non-boundary species (147 - 3 = 144).
    assert m.getNumRules() == 144


def test_simulate_via_tellurium_normal_regime(tmp_path):
    """End-to-end: export → load via tellurium → simulate Normal regime."""
    te = pytest.importorskip("tellurium")
    out = tmp_path / "model.xml"
    net = load_network(DATA)
    export_sbml(net, out)

    r = te.loadSBMLModel(str(out))
    r["Hypo"] = 0.01
    r["NL"]   = 0.80
    r["HL"]   = 0.01
    res = r.simulate(0, 100, 2)

    # Under Normal regime: SOX9 should be high, ROS low.
    sox9_final = res["[SOX9]"][-1]
    ros_final = res["[ROS]"][-1]
    assert sox9_final > 0.5, f"SOX9 should be high under Normal, got {sox9_final}"
    assert ros_final < 0.1,  f"ROS should be low under Normal, got {ros_final}"
