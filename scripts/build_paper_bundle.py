"""Bundle the manuscript's results for the webapp's Paper results tab.

Reads the CSV tables written by the run_*.py scripts and writes
results/web_bundle/paper_results.json: one section per manuscript figure
(topology, baseline panels, transition heatmaps, falsification, rescue
rankings and node-resolved rescue responses).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

from np_mt_rnm.categories import CATEGORY_LABELS, nodes_in_category
from np_mt_rnm.transitions import HYPO_TO_NORMAL_PATH, NORMAL_TO_HYPER_PATH

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
OUT = ROOT / "results" / "web_bundle" / "paper_results.json"

# Manuscript figure panels, in the manuscript's order and wording.
BASELINE_FIGURES = [
    {
        "id": "BSL_ECM_GF",
        "caption": "Baseline steady-state responses under Hypo, Normal and Hyper loading (manuscript Fig. 4).",
        "panels": [
            ("ecm_matrix", "ECM anabolism & phenotype markers"),
            ("growth_factor", "Growth factors"),
            ("transcription_factor", "Transcription Factors"),
            ("cytokines_chemokines_proteases", "Cytokines, chemokines, proteases & others"),
        ],
    },
    {
        "id": "BSL_MET_OXI",
        "caption": "Baseline steady-state responses under Hypo, Normal and Hyper loading (manuscript Fig. 5).",
        "panels": [
            ("metabolic", "Metabolic & related"),
            ("ion_channel", "Ion channels & related"),
            ("oxidative_proteostasis", "Oxidative-stress defense & proteostasis"),
            ("cell_fate", "Cell survival, apoptosis & mitophagy/DNA-damage"),
        ],
    },
]

# The six functional modules of the transition and rescue analyses.
MODULES = [
    ("ecm_matrix", "ECM anabolism & phenotype markers", "ECM"),
    ("growth_factor", "Growth factors", "GF"),
    ("transcription_factor", "Transcription factors", "TF"),
    ("cytokines_chemokines_proteases", "Cytokines, chemokines, proteases & related mediators", "CYT"),
    ("oxidative_proteostasis", "Oxidative-stress defense & proteostasis", "OXI"),
    ("cell_fate", "Cell survival, apoptosis, mitophagy & DNA damage", "CSF"),
]

REPRESENTATIVE_NODES = [
    "COL2A1", "ACAN", "TIMP3", "COL1A1", "COL10A1", "SOX9",
    "NF-κB", "ROS", "IGF1", "VEGF", "Bcl2", "CASP3",
]
REPRESENTATIVE_ANABOLIC = ["COL2A1", "ACAN", "TIMP3", "SOX9", "IGF1", "Bcl2"]


def _clean(x):
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    return x


def main() -> None:
    baseline = pd.read_csv(TABLES / "baseline_summary.csv")
    means = baseline.pivot(index="node", columns="regime", values="mean_activation")
    node_set = set(means.index)

    def members(cat):
        return [n for n in nodes_in_category(cat) if n in node_set]

    topo = pd.read_csv(TABLES / "topology_metrics.csv")
    topology = [
        {
            "node": r.node,
            "out_activating": int(r.out_activating),
            "out_inhibiting": int(r.out_inhibiting),
            "betweenness": float(r.betweenness),
            "harmonic_closeness": float(r.harmonic_closeness),
        }
        for r in topo.itertuples()
    ]

    baseline_figs = []
    for fig in BASELINE_FIGURES:
        panels = []
        for cat, title in fig["panels"]:
            # The manuscript's metabolic panel follows NP_MT_RNM_FALSIFY4_1.m,
            # whose group list has no NutD.
            nodes = [n for n in members(cat) if not (cat == "metabolic" and n == "NutD")]
            panels.append({
                "category": cat,
                "title": title,
                "nodes": nodes,
                "series": {
                    reg: [float(means.loc[n, reg]) for n in nodes]
                    for reg in ("Hypo", "Normal", "Hyper")
                },
            })
        baseline_figs.append({"id": fig["id"], "caption": fig["caption"], "panels": panels})

    transitions = {}
    for key, path in [("hypo_to_normal", HYPO_TO_NORMAL_PATH), ("normal_to_hyper", NORMAL_TO_HYPER_PATH)]:
        df = pd.read_csv(TABLES / f"transition_{key}.csv")
        transitions[key] = {
            "steps": [dict(s) for s in path],
            "means": {n: [float(v) for v in df[n]] for n in df.columns if n != "step"},
        }

    fals = pd.read_csv(TABLES / "falsification_rules.csv")
    falsification = [
        {
            "node": r.node,
            "class": r._2,
            "delta": float(r.delta),
            "ci_lower": float(r.ci_lower),
            "ci_upper": float(r.ci_upper),
            "passed": bool(r.passed),
        }
        for r in fals.itertuples()
    ]

    rankings = pd.read_csv(TABLES / "true_rescue_rankings.csv")
    per_node = pd.read_csv(TABLES / "rescue_per_node_deltas.csv").set_index("perturbation")
    rescue = []
    for cat, title, tag in MODULES:
        rk = rankings[rankings.category == cat].sort_values("rank")
        nodes = [n for n in members(cat) if n in per_node.columns]
        rescue.append({
            "category": cat,
            "title": title,
            "tag": tag,
            "top20": [
                {"strategy": r.rescue_strategy, "rescue_percent": float(r.rescue_percent)}
                for r in rk.itertuples()
            ],
            "nodes": nodes,
            "perturbations": list(per_node.index),
            # mean paired Δ (perturbed − Hyper), rows = perturbations
            "mean_delta": [[float(per_node.loc[p, n]) for n in nodes] for p in per_node.index],
        })

    bundle = {
        "provenance": {
            "source": "np-mt-rnm results/tables (scripts/run_all.py)",
            "network": "MT_PRIMARY4_1.xlsx, 147 nodes, 357 edges",
            "seeding": "MATLAB twister seeds of RESCUE_NEW4_1_final.m (Hypo 1001+r, Normal 2001+r, Hyper 3001+r)",
        },
        "topology": topology,
        "baseline_figures": baseline_figs,
        "transitions": transitions,
        "transition_modules": [
            {"category": c, "title": t, "tag": g, "nodes": members(c)} for c, t, g in MODULES
        ],
        "representative": {"nodes": REPRESENTATIVE_NODES, "anabolic": REPRESENTATIVE_ANABOLIC},
        "falsification": {"ftol": 0.02, "n_boot": 10000, "rules": falsification},
        "rescue": rescue,
        "category_labels": CATEGORY_LABELS,
    }
    OUT.write_text(json.dumps(bundle, ensure_ascii=False, default=_clean))
    print(f"[paper-bundle] wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
