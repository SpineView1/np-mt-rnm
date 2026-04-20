"""Convert results/ CSVs + arrays into webapp-consumable JSON.

Outputs to results/web_bundle/*.json. The webapp (np-mt-rnm-web) reads
these JSON files directly — it never parses the CSVs.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from np_mt_rnm.categories import CATEGORY_LABELS, get_primary_category
from np_mt_rnm.network import load_network
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

    (BUNDLE / "network.json").write_text(json.dumps(_build_network_json(net)))
    (BUNDLE / "baseline.json").write_text(json.dumps(_build_baseline_json()))
    (BUNDLE / "falsification.json").write_text(json.dumps(_build_falsification_json()))
    (BUNDLE / "rescue.json").write_text(json.dumps(_build_rescue_json()))
    (BUNDLE / "transitions.json").write_text(json.dumps(_build_transitions_json()))

    sizes = {p.name: p.stat().st_size for p in sorted(BUNDLE.glob("*.json"))}
    print(f"[web_bundle] wrote 5 JSON files to {BUNDLE}")
    for name, size in sizes.items():
        print(f"  {name:25s} {size/1024:.1f} KB")


if __name__ == "__main__":
    main()
