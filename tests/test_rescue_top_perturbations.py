"""The committed rescue tables carry the manuscript's top intervention per module.

Supplementary Table "Highest-scoring intervention for different functional modules".
"""
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

TOP = {
    "ecm_matrix": ("SOX9-ROS", 85.9),
    "oxidative_proteostasis": ("FOXO-ROS", 68.9),
    "growth_factor": ("SOX9-NF-κB", 66.0),
    "cytokines_chemokines_proteases": ("SOX9-NF-κB", 57.0),
    "cell_fate": ("NRF2-PIEZO1", 37.5),
    "transcription_factor": ("SOX9-ROS", 30.7),
    "metabolic": ("PPARγ-ROS", 23.3),
    "ion_channel": ("PPARγ-PIEZO1", 10.4),
}


@pytest.mark.parametrize("cat", list(TOP))
def test_module_top_matches_manuscript(cat):
    df = pd.read_csv(ROOT / "results" / "tables" / "rescue_module_top.csv").set_index("category")
    label, pct = TOP[cat]
    assert df.loc[cat, "top_intervention"] == label
    assert f"{df.loc[cat, 'rescue_percent']:.1f}" == f"{pct:.1f}"


def test_phenotype_panel_table_matches_manuscript():
    df = pd.read_csv(ROOT / "results" / "tables" / "rescue_phenotype_panels.csv")
    df = df.set_index("regime_or_intervention")
    expected = {
        "NL baseline": (0.985, 0.001), "HL baseline": (0.059, 0.973),
        "SOX9-ROS": (0.589, 0.653), "SOX9-NF-κB": (0.718, 0.371),
        "SOX9-P38": (0.561, 0.654), "SOX9-mTORC1": (0.562, 0.655),
        "SOX9-PIEZO1": (0.557, 0.645), "SOX9-RhoA-E": (0.556, 0.628),
        "SOX9": (0.556, 0.656), "PPARγ-ROS": (0.463, 0.655),
        "PPARγ-P38": (0.434, 0.658), "PPARγ-PIEZO1": (0.431, 0.649),
        "PPARγ-NF-κB": (0.458, 0.445), "PPARγ-RhoA-E": (0.430, 0.634),
    }
    got = {k: (round(df.loc[k, "anabolic_mean"], 3), round(df.loc[k, "catabolic_mean"], 3))
           for k in expected}
    assert got == expected
