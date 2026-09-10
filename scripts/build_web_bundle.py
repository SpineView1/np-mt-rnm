"""Convert results/ CSVs + arrays into webapp-consumable JSON.

Outputs to results/web_bundle/*.json. The webapp (np-mt-rnm-web) reads
these JSON files directly — it never parses the CSVs.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from np_mt_rnm.categories import (
    CATEGORY_LABELS,
    CATEGORY_ORDER,
    get_primary_category,
    nodes_in_category,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.figures import REGIME_COLORS
from np_mt_rnm.simulation import PAPER_REGIMES
from np_mt_rnm.topology import compute_topology

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "MT_PRIMARY4_1.xlsx"
RESULTS = ROOT / "results"
BUNDLE = RESULTS / "web_bundle"


def _build_network_json(net) -> dict:
    metrics = compute_topology(net)
    nodes = []
    for name in net.node_names:
        try:
            primary = get_primary_category(name)
        except KeyError:
            primary = "other"
        nodes.append({
            "id": name,
            "category": primary,
            "category_label": CATEGORY_LABELS.get(primary, primary),
            "out_degree": int(metrics.signed_out_degree.get(name, 0)),
            "betweenness": float(metrics.betweenness.get(name, 0.0)),
        })
    edges = []
    for i in range(len(net.node_names)):
        for j in range(len(net.node_names)):
            if net.mact[i, j]:
                edges.append({
                    "source": net.node_names[j],
                    "target": net.node_names[i],
                    "sign": 1,
                })
            if net.minh[i, j]:
                edges.append({
                    "source": net.node_names[j],
                    "target": net.node_names[i],
                    "sign": -1,
                })
    return {"nodes": nodes, "edges": edges}


def _build_baseline_json() -> dict:
    df = pd.read_csv(RESULTS / "tables" / "baseline_summary.csv")
    out: dict[str, dict[str, dict[str, float]]] = {}
    for regime in df["regime"].unique():
        sub = df[df["regime"] == regime]
        out[regime] = {
            row["node"]: {
                "mean": float(row["mean_activation"]),
                "std": float(row["std_activation"]),
            }
            for _, row in sub.iterrows()
        }
    return out


# The four panels of paper Figure 4, in caption order:
#   (A) ECM anabolism and phenotype markers
#   (B) Growth factors
#   (C) Transcription factors
#   (D) Cytokines, chemokines, proteases, and related mediators
FIGURE4_PANELS: tuple[str, ...] = (
    "ecm_matrix",
    "growth_factor",
    "transcription_factor",
    "cytokines_chemokines_proteases",
)

# Panel headings worded as the paper's Figure 4 caption, which differs slightly
# from the internal CATEGORY_LABELS used elsewhere in the app.
FIGURE4_TITLES: dict[str, str] = {
    "ecm_matrix": "ECM anabolism and phenotype markers",
    "growth_factor": "Growth factors",
    "transcription_factor": "Transcription factors",
    "cytokines_chemokines_proteases":
        "Cytokines, chemokines, proteases, and related mediators",
}


def _build_baseline_groups_json(net, baseline: dict) -> dict:
    """Pre-shape the baseline means into per-category panels for the webapp.

    The webapp plots grouped bars per biological category, one series per
    loading regime, exactly as plot_grouped_bars_no_errors does in
    NP_MT_RNM_FALSIFY4_1.m. Doing the grouping here keeps the node ordering
    and category membership in one place instead of duplicating the MATLAB
    group lists in JavaScript.
    """
    present = set(net.node_names)
    groups: dict[str, dict] = {}

    for cat in CATEGORY_ORDER:
        if cat == "other":
            continue
        nodes = [n for n in nodes_in_category(cat) if n in present]
        if not nodes:
            continue
        series = {}
        for regime in PAPER_REGIMES:
            per_node = baseline[regime]
            series[regime] = {
                "mean": [per_node[n]["mean"] for n in nodes],
                "std": [per_node[n]["std"] for n in nodes],
            }
        groups[cat] = {
            "label": CATEGORY_LABELS[cat],
            "nodes": nodes,
            "series": series,
        }

    missing_panels = [c for c in FIGURE4_PANELS if c not in groups]
    if missing_panels:
        raise ValueError(f"Figure 4 panels missing from bundle: {missing_panels}")

    return {
        "regimes": list(PAPER_REGIMES),
        "regime_colors": {r: REGIME_COLORS[r] for r in PAPER_REGIMES},
        "figure4_panels": list(FIGURE4_PANELS),
        "figure4_titles": dict(FIGURE4_TITLES),
        "groups": groups,
        "provenance": {
            "source": "np-mt-rnm results/tables/baseline_summary.csv",
            "statistic": "mean activation over 100 replicates from random initial conditions",
            "note": (
                "Means only, matching the paper's Figure 4 "
                "(plot_grouped_bars_no_errors draws no error bars). Standard "
                "deviations are included per node because several nodes are "
                "bimodal across replicates under Hyper loading, where a mean "
                "describes no individual steady state."
            ),
        },
    }


def _build_falsification_json() -> dict:
    df = pd.read_csv(RESULTS / "tables" / "falsification_rules.csv")
    return {"rules": df.to_dict(orient="records")}


def _build_rescue_json() -> dict:
    summary = pd.read_csv(RESULTS / "tables" / "rescue_summary.csv")
    per_node = pd.read_csv(RESULTS / "tables" / "rescue_per_node_deltas.csv")
    return {
        "summary": summary.to_dict(orient="records"),
        "per_node": per_node.to_dict(orient="records"),
    }


def _build_transitions_json() -> dict:
    out = {}
    for name in ("hypo_to_normal", "normal_to_hyper"):
        df = pd.read_csv(RESULTS / "tables" / f"transition_{name}.csv")
        out[name] = df.to_dict(orient="records")
    return out


def main() -> None:
    BUNDLE.mkdir(parents=True, exist_ok=True)
    net = load_network(DATA)

    baseline = _build_baseline_json()
    (BUNDLE / "network.json").write_text(json.dumps(_build_network_json(net)))
    (BUNDLE / "baseline.json").write_text(json.dumps(baseline))
    (BUNDLE / "baseline_groups.json").write_text(
        json.dumps(_build_baseline_groups_json(net, baseline))
    )
    (BUNDLE / "falsification.json").write_text(json.dumps(_build_falsification_json()))
    (BUNDLE / "rescue.json").write_text(json.dumps(_build_rescue_json()))
    (BUNDLE / "transitions.json").write_text(json.dumps(_build_transitions_json()))

    sizes = {p.name: p.stat().st_size for p in sorted(BUNDLE.glob("*.json"))}
    print(f"[web_bundle] wrote {len(sizes)} JSON files to {BUNDLE}")
    for name, size in sizes.items():
        print(f"  {name:25s} {size/1024:.1f} KB")


if __name__ == "__main__":
    main()
