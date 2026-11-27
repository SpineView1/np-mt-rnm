"""The rescue screen reproduces the manuscript to the printed decimal.

Every value below is printed in the manuscript (Workineh, Chemorion & Noailly,
Results, "Hyper-to-Normal rescue screening via node clamping", and the
Supplementary rescue tables), which RESCUE_NEW4_1_48_targets_HQ.m writes with
sprintf('%.1f%%') / '%.3f'.

The network is the 353-interaction, 146-node corpus (Supplementary Data 1).
Initial states are MATLAB's twister draws with the rescue script's
per-replicate seeds, so the same numbers come out of the MATLAB script.
"""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import (
    PAPER_RESCUE_MODULES,
    RESCUE_GROUPS,
    enumerate_perturbations,
    rank_group_rescue,
    run_perturbation,
)
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates

ROOT = Path(__file__).resolve().parents[1]

# module -> {manuscript label: printed rescue %}
PAPER = {
    "ecm_matrix": {
        "SOX9-ROS": 85.9, "SOX9-NF-κB": 85.2, "SOX9-P38": 84.2, "SOX9-mTORC1": 83.8,
        "SOX9-PIEZO1": 83.0, "SOX9-RhoA-E": 82.9, "SOX9": 82.8, "PPARγ": 48.8,
        "PPARγ-ROS": 52.3,
    },
    "growth_factor": {
        "SOX9-NF-κB": 66.0, "PPARγ-NF-κB": 24.1, "SOX9-RhoA-E": 22.5, "SOX9": 20.7,
    },
    "transcription_factor": {
        "SOX9-ROS": 30.7, "SOX9-NF-κB": 23.8, "PPARγ-ROS": 23.4, "SOX9": 15.4, "ROS": 12.6,
    },
    "cytokines_chemokines_proteases": {
        "SOX9-NF-κB": 57.0, "PPARγ-NF-κB": 46.2, "SOX9-RhoA-E": 20.8,
        "PPARγ-RhoA-E": 20.2, "SOX9": 17.9, "PPARγ": 17.6, "NF-κB": 3.6,
    },
    "metabolic": {"PPARγ-ROS": 23.3, "FOXO-ROS": 17.8, "SIRT1-ROS": 15.9, "PPARγ": 14.6},
    "ion_channel": {
        "PPARγ-PIEZO1": 10.4, "SIRT1-PIEZO1": 10.4, "FOXO-PIEZO1": 10.4, "PIEZO1": 6.0,
        "PPARγ-ROS": 8.7, "SIRT1-ROS": 8.7, "FOXO-ROS": 8.7,
    },
    "oxidative_proteostasis": {
        "FOXO-ROS": 68.9, "SIRT1-ROS": 60.8, "NRF2": 60.4, "FOXO": 24.1, "ROS": 27.8,
    },
    "cell_fate": {
        "NRF2-PIEZO1": 37.5, "NRF2": 31.0, "PPARγ-ROS": 34.9, "FOXO-ROS": 34.7, "ROS": 34.2,
    },
}

# Highest-ranked strategy per module, as stated in the manuscript text.
PAPER_TOP = {
    "ecm_matrix": "SOX9-ROS",
    "growth_factor": "SOX9-NF-κB",
    "transcription_factor": "SOX9-ROS",
    "cytokines_chemokines_proteases": "SOX9-NF-κB",
    "metabolic": "PPARγ-ROS",
    "oxidative_proteostasis": "FOXO-ROS",
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
    names = list(net.node_names)
    final = np.column_stack([r.perturbed_mean for r in results])
    labels = [r.perturbation.short_label for r in results]
    rankings = {
        cat: rank_group_rescue(
            group=cat, group_nodes=list(RESCUE_GROUPS[cat]), node_names=names,
            hyper_mean=hyper.mean(), normal_mean=normal.mean(),
            final_states=final, strategy_labels=labels, top_n=48,
        )
        for cat in PAPER_RESCUE_MODULES
    }
    return {"rankings": rankings, "names": names, "final": final, "labels": labels,
            "hyper": hyper.mean()}


@pytest.mark.parametrize("cat", list(PAPER))
def test_rescue_percentages_match_manuscript(screen, cat):
    ranking = screen["rankings"][cat]
    got = dict(zip(ranking.strategy_labels, ranking.rescue_percent))
    mismatches = {
        lab: (printed, round(got[lab], 2))
        for lab, printed in PAPER[cat].items()
        if f"{got[lab]:.1f}" != f"{printed:.1f}"
    }
    assert not mismatches, mismatches


@pytest.mark.parametrize("cat", list(PAPER_TOP))
def test_top_strategy_matches_manuscript(screen, cat):
    assert screen["rankings"][cat].strategy_labels[0] == PAPER_TOP[cat]


def test_ecm_node_resolved_values(screen):
    """SOX9-ROS raises ACAN, COL2A1, TIMP3 from 0.216, 0.252, <0.001 to 0.986, 0.994, 0.716."""
    names, final, labels, hyper = (screen[k] for k in ("names", "final", "labels", "hyper"))
    j = labels.index("SOX9-ROS")
    ix = [names.index(n) for n in ("ACAN", "COL2A1", "TIMP3")]
    assert [round(v, 3) for v in hyper[ix[:2]]] == [0.216, 0.252]
    assert hyper[ix[2]] < 0.001
    assert [round(v, 3) for v in final[ix, j]] == [0.986, 0.994, 0.716]
