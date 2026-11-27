"""Bundle the manuscript's results for the webapp's Paper results tab.

Reads the CSV tables written by the run_*.py scripts and writes
results/web_bundle/paper_results.json: one section per manuscript figure
(topology, baseline panels, falsification, rescue rankings and node-resolved
rescue responses), plus the falsification benchmark registry and the rescue
phenotype-panel table.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

from np_mt_rnm import paper
from np_mt_rnm.categories import CATEGORY_LABELS, nodes_in_category
from np_mt_rnm.falsification import FALS_ALPHA, FALS_NBOOT, outcome_counts  # noqa: F401
from np_mt_rnm.rescue import (
    ANABOLIC_PANEL,
    ANABOLIC_UP_NODES,
    CATABOLIC_DOWN_NODES,
    CATABOLIC_PANEL,
    RESCUE_GROUPS,
)

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
OUT = ROOT / "results" / "web_bundle" / "paper_results.json"

# Manuscript figure panels, in the manuscript's order and wording.
BASELINE_FIGURES = [
    {
        "id": "BSL_ECM_GF",
        "panels": [
            ("ecm_matrix", "ECM anabolism & phenotype markers"),
            ("growth_factor", "Growth factors"),
            ("transcription_factor", "Transcription Factors"),
            ("cytokines_chemokines_proteases", "Cytokines, chemokines, proteases & others"),
        ],
    },
    {
        "id": "BSL_MET_OXI",
        "panels": [
            ("metabolic", "Metabolic & related"),
            ("ion_channel", "Ion channels & related"),
            ("oxidative_proteostasis", "Oxidative-stress defense & proteostasis"),
            ("cell_fate", "Cell survival, apoptosis & mitophagy/DNA-damage"),
        ],
    },
]

# The eight rescue modules of manuscript Figures 6 (A-D) and 7 (A-D).
MODULES = [
    ("ecm_matrix", "ECM anabolism & phenotype markers", "ECM", "Figure 6", "A"),
    ("growth_factor", "Growth factors", "GF", "Figure 6", "B"),
    ("transcription_factor", "Transcription factors", "TF", "Figure 6", "C"),
    ("cytokines_chemokines_proteases", "Cytokines, chemokines, proteases & related mediators",
     "CYT", "Figure 6", "D"),
    ("metabolic", "Metabolic & related", "MET", "Figure 7", "A"),
    ("ion_channel", "Ion channels & related", "ION", "Figure 7", "B"),
    ("oxidative_proteostasis", "Oxidative-stress defense & proteostasis", "OXI", "Figure 7", "C"),
    ("cell_fate", "Cell survival, apoptosis, mitophagy & DNA damage", "CS", "Figure 7", "D"),
]


def _clean(x):
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    return x


def main() -> None:
    ci = pd.read_csv(TABLES / "baseline_means_95ci.csv").set_index("Node")
    node_set = set(ci.index)
    col = {"Hypo": "Hypo", "Normal": "NL", "Hyper": "HL"}

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

    stats = pd.read_csv(TABLES / "baseline_nonparametric_statistics.csv").set_index("Node")
    baseline_figs = []
    for fig in BASELINE_FIGURES:
        panels = []
        for cat, title in fig["panels"]:
            nodes = members(cat)
            panel = {"category": cat, "title": title, "nodes": nodes}
            for key, suffix in (("series", "Mean"), ("sd", "SD"),
                                ("ci_low", "CI_Low"), ("ci_high", "CI_High")):
                panel[key] = {
                    reg: [float(ci.loc[n, f"{c}_{suffix}"]) for n in nodes]
                    for reg, c in col.items()
                }
            panel["q_nl_hl"] = [_clean(float(stats.loc[n, "q_NL_HL"])) for n in nodes]
            panel["q_kw"] = [_clean(float(stats.loc[n, "q_KW_3group"])) for n in nodes]
            panels.append(panel)
        baseline_figs.append({"id": fig["id"], "panels": panels})

    fals = pd.read_csv(TABLES / "falsification_rules.csv", keep_default_na=False)
    rules = [
        {
            "node": r["node"],
            "class": r["class"],
            "direction": r["direction"],
            "literature_direction": r["literature_direction"],
            "benchmark_ids": r["benchmark_ids"].split(";"),
            "references": r["reference"].split(";"),
            "mean_normal": float(r["mean_normal"]),
            "mean_hyper": float(r["mean_hyper"]),
            "delta": float(r["delta"]),
            "ci_lower": float(r["ci_lower"]),
            "ci_upper": float(r["ci_upper"]),
            "outcome": r["outcome"],
            "passed": str(r["passed"]) == "True",
        }
        for _, r in fals.iterrows()
    ]
    counts = {k: int((fals["outcome"] == k).sum())
              for k in ("Concordant", "Contradictory", "Indeterminate")}
    # Baseline CIs of the benchmark nodes for the NL / HL profile panel (Fig. 5A).
    for r in rules:
        r["nl_ci"] = [float(ci.loc[r["node"], "NL_CI_Low"]), float(ci.loc[r["node"], "NL_CI_High"])]
        r["hl_ci"] = [float(ci.loc[r["node"], "HL_CI_Low"]), float(ci.loc[r["node"], "HL_CI_High"])]
        r["nl_sd"] = float(ci.loc[r["node"], "NL_SD"])
        r["hl_sd"] = float(ci.loc[r["node"], "HL_SD"])

    rankings = pd.read_csv(TABLES / "true_rescue_rankings.csv")
    per_node = pd.read_csv(TABLES / "rescue_per_node_deltas.csv").set_index("perturbation")
    per_sd = pd.read_csv(TABLES / "rescue_per_node_sd.csv").set_index("perturbation")
    rescue = []
    for cat, title, tag, fig_no, panel in MODULES:
        rk = rankings[rankings.category == cat].sort_values("rank")
        nodes = [n for n in RESCUE_GROUPS[cat] if n in per_node.columns]
        top2 = list(rk["rescue_strategy"].head(2))
        rescue.append({
            "category": cat,
            "title": title,
            "tag": tag,
            "figure": fig_no,
            "panel": panel,
            "figure_file": f"Top20_{tag}_TrueRescue",
            "top20": [
                {"strategy": r.rescue_strategy, "rescue_percent": float(r.rescue_percent)}
                for r in rk.head(20).itertuples()
            ],
            "all48": {r.rescue_strategy: float(r.rescue_percent) for r in rk.itertuples()},
            "nodes": nodes,
            "perturbations": list(per_node.index),
            # mean / SD of the paired Δ (perturbed − Hyper), rows = perturbations
            "mean_delta": [[float(per_node.loc[p, n]) for n in nodes] for p in per_node.index],
            "sd_delta": [[float(per_sd.loc[p, n]) for n in nodes] for p in per_node.index],
            "top2": top2,
        })

    pheno = pd.read_csv(TABLES / "rescue_phenotype_panels.csv")
    module_top = pd.read_csv(TABLES / "rescue_module_top.csv")

    bundle = {
        "provenance": {
            "source": "np-mt-rnm results/tables (scripts/run_all.py)",
            "network": "MT_PRIMARY4_1.xlsx, 146 nodes, 353 interactions",
            "seeding": "MATLAB twister seeds of RESCUE_NEW4_1_48_targets_HQ.m "
                       "(Hypo 1001+r, Normal 2001+r, Hyper 3001+r)",
            "matlab": ["legacy/NP_MT_FALSIFY4_1_39.m", "legacy/RESCUE_NEW4_1_48_targets_HQ.m"],
        },
        "topology": topology,
        "baseline_figures": baseline_figs,
        "baseline_note": (
            "Error bars: 95% percentile bootstrap CIs of the mean (10,000 resamples), "
            "as drawn by NP_MT_FALSIFY4_1_39.m. Per-node SDs are included as well. "
            "q-values: Kruskal–Wallis and NL–HL Wilcoxon rank-sum, BH-FDR across nodes."
        ),
        "falsification": {
            "criterion": "zero-based 95% bootstrap CI",
            "n_boot": FALS_NBOOT,
            "alpha": FALS_ALPHA,
            "counts": counts,
            "rules": rules,
        },
        "rescue_design": {
            "suppression": list(CATABOLIC_DOWN_NODES),
            "activation": list(ANABOLIC_UP_NODES),
            "n_singles": len(CATABOLIC_DOWN_NODES) + len(ANABOLIC_UP_NODES),
            "n_duals": len(CATABOLIC_DOWN_NODES) * len(ANABOLIC_UP_NODES),
            "anabolic_panel": list(ANABOLIC_PANEL),
            "catabolic_panel": list(CATABOLIC_PANEL),
        },
        "rescue": rescue,
        "rescue_module_top": module_top.to_dict(orient="records"),
        "phenotype_panels": pheno.to_dict(orient="records"),
        "category_labels": CATEGORY_LABELS,
        # Manuscript figure numbers and captions, verbatim (np_mt_rnm/paper.py).
        "figures": {k: {"number": n, "caption": c} for k, (n, c) in paper.FIGURES.items()},
        "node_resolved_intro": paper.NODE_RESOLVED_INTRO,
        "falsification_summary": paper.FALSIFICATION_SUMMARY,
        "rescue_summary": paper.RESCUE_SUMMARY,
        "paper": {
            "title": paper.TITLE,
            "authors": list(paper.AUTHORS),
            "affiliation": paper.AFFILIATION,
            "keywords": list(paper.KEYWORDS),
            "abstract": paper.ABSTRACT,
        },
    }
    OUT.write_text(json.dumps(bundle, ensure_ascii=False, default=_clean))
    print(f"[paper-bundle] wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
