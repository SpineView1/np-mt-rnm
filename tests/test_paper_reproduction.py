"""The rescue screen reproduces the manuscript to the printed decimal.

Every value below is printed in the manuscript (Workineh, Chemorion & Noailly,
Section 3.5) or as a bar label in its Top-20 rescue figures (ECM_rescue.png,
GF_rescue.png, TF_rescue.png, CYT_rescue.png, OXI_rescue.png, CSF_rescue.png),
which RESCUE_NEW4_1_final.m writes with sprintf('%.1f%%').

The match is exact only when (i) initial states are MATLAB's twister draws with
the rescue script's per-replicate seeds and (ii) the network is the 357-edge
MT_PRIMARY4_1.xlsx. On the 353-edge 2026-09 revision the TF and OXI modules
move by 0.1-0.3 points.
"""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.categories import nodes_in_category
from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import enumerate_perturbations, rank_group_rescue, run_perturbation
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates

ROOT = Path(__file__).resolve().parents[1]


def _label(s: str) -> str:
    """'SOX9-FAK-E' (figure label) -> 'SOX9↑ + FAK-E↓' (port label)."""
    kd = ("RhoA-E", "PIEZO1", "PI3K-E", "FAK-E", "ROS")
    if s in kd:
        return f"{s}↓"
    for k in kd:
        if s.endswith("-" + k):
            return f"{s[: -len(k) - 1]}↑ + {k}↓"
    return f"{s}↑"


# category id -> {figure label: printed rescue %}
PAPER = {
    "ecm_matrix": {"SOX9-FAK-E": 90.5, "SOX9-ROS": 85.9, "SOX9": 82.8},
    "growth_factor": {"SOX9-FAK-E": 31.5, "SOX9": 20.7},
    "transcription_factor": {
        "SOX9-ROS": 30.6, "PPARγ-ROS": 23.2, "SOX9-PIEZO1": 18.4, "IκBα-ROS": 17.4,
        "SOX9-PI3K-E": 15.8, "SOX9": 15.5, "SOX9-FAK-E": 15.0, "NRF2-PIEZO1": 14.2,
        "ROS": 12.4, "SOX9-RhoA-E": 12.4, "HIF-1α-ROS": 12.3, "PPARγ-PIEZO1": 12.0,
        "NRF2-PI3K-E": 11.7, "NRF2-FAK-E": 11.4, "NRF2": 11.3, "NRF2-ROS": 11.3,
        "PPARγ-FAK-E": 10.0, "PPARγ-PI3K-E": 9.3, "PPARγ": 9.3, "NRF2-RhoA-E": 8.4,
    },
    "cytokines_chemokines_proteases": {"SOX9-FAK-E": 22.4},
    "oxidative_proteostasis": {
        "NRF2-PIEZO1": 60.3, "NRF2-RhoA-E": 60.3, "NRF2": 60.3, "NRF2-ROS": 60.3,
        "NRF2-PI3K-E": 60.3, "NRF2-FAK-E": 59.0, "HIF-1α-FAK-E": 34.5,
        "PPARγ-ROS": 28.1, "HIF-1α-ROS": 27.7, "SOX9-ROS": 27.7, "IκBα-ROS": 27.7,
        "ROS": 27.7, "FAK-E": 23.7, "IκBα-FAK-E": 23.7, "SOX9-FAK-E": 23.7,
        "PPARγ-FAK-E": 23.7, "HIF-1α-RhoA-E": 8.4, "HIF-1α-PI3K-E": 8.3,
        "HIF-1α": 8.0, "HIF-1α-PIEZO1": 8.0,
    },
    "cell_fate": {"NRF2-PIEZO1": 37.5, "PPARγ-ROS": 34.9, "ROS": 34.2},
}

# Highest-ranked strategy per module, as stated in the manuscript text.
PAPER_TOP = {
    "ecm_matrix": "SOX9-FAK-E",
    "growth_factor": "SOX9-FAK-E",
    "transcription_factor": "SOX9-ROS",
    "cytokines_chemokines_proteases": "SOX9-FAK-E",
    "cell_fate": "NRF2-PIEZO1",
}


@pytest.fixture(scope="module")
def screen():
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    hyper = run_replicates(
        net, REGIME_PRESETS["Hyper"], n_reps=100, seed=REGIME_SEEDS["Hyper"], n_jobs=-1
    )
    normal = run_replicates(
        net, REGIME_PRESETS["Normal"], n_reps=100, seed=REGIME_SEEDS["Normal"], n_jobs=-1
    )
    results = [
        run_perturbation(
            net, p.anabolic_up, p.catabolic_down, n_reps=100, n_jobs=-1, baseline=hyper
        )
        for p in enumerate_perturbations()
    ]
    rankings = {}
    for cat in PAPER:
        rankings[cat] = rank_group_rescue(
            group=cat,
            group_nodes=nodes_in_category(cat),
            node_names=list(net.node_names),
            hyper_mean=hyper.mean(),
            normal_mean=normal.mean(),
            final_states=np.column_stack([r.perturbed_mean for r in results]),
            strategy_labels=[r.perturbation.label for r in results],
            top_n=35,
        )
    return rankings


@pytest.mark.parametrize("cat", list(PAPER))
def test_rescue_percentages_match_manuscript(screen, cat):
    ranking = screen[cat]
    got = dict(zip(ranking.strategy_labels, ranking.rescue_percent))
    mismatches = {
        fig: (printed, round(got[_label(fig)], 2))
        for fig, printed in PAPER[cat].items()
        if f"{got[_label(fig)]:.1f}" != f"{printed:.1f}"
    }
    assert not mismatches, mismatches


@pytest.mark.parametrize("cat", list(PAPER_TOP))
def test_top_strategy_matches_manuscript(screen, cat):
    assert screen[cat].strategy_labels[0] == _label(PAPER_TOP[cat])
