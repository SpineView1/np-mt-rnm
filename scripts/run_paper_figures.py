"""Draw every computational figure of the manuscript, named as in the paper.

Reads the outputs of run_topology / run_baseline / run_falsification /
run_rescue (no re-simulation) and writes results/figures/paper/<Name>.png:

  TOPO_STATS.png                      Figure 2
  BSL_ECM_GF.png, BSL_MET_OXI.png     Figures 3-4 (95% bootstrap CIs of the mean)
  NP_MT_FALS_RECT.png, NP_MT_FALS.png Figure 5 (A-B, C)
  Top20_<tag>_TrueRescue.png          Figures 6-7 (8 modules)
  NodeResolved_ECM_GF_TF_CYT.png,
  NodeResolved_MET_ION_SURV_OXI.png   Supplementary node-resolved figures
  heatmaps/heatmap_{all,singles,duals}_<module>.png  rescue Δ with BH-FDR stars
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm import paper_figures as pf
from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import N_SINGLES, RESCUE_GROUPS
from np_mt_rnm.statistics import bh_fdr, permutation_pvalues

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
REPS = RESULTS / "replicates"
OUT = RESULTS / "figures" / "paper"

# rng(PERMUTATION_SEED, 'twister') in RESCUE_NEW4_1_48_targets_HQ.m, BASE_SEED + 4000.
PERMUTATION_SEED = 4001
N_PERM = 1000


def topology() -> None:
    df = pd.read_csv(TABLES / "topology_metrics.csv")
    n = len(df)
    # The manuscript plots raw values: betweenness as a path count and
    # harmonic closeness as the sum of 1/d(i,j) over outgoing shortest paths.
    pf.plot_topo_stats(
        df["node"], df["out_activating"], df["out_inhibiting"],
        df["betweenness"] * (n - 1) * (n - 2),
        df["harmonic_closeness"] * (n - 1),
        OUT / "TOPO_STATS.png",
        node_order=pf.edge_list_node_order(load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")),
    )


def _baseline_ci():
    t = pd.read_csv(TABLES / "baseline_means_95ci.csv")
    names = list(t["Node"])
    col = {"Hypo": "Hypo", "Normal": "NL", "Hyper": "HL"}
    means = {r: t[f"{c}_Mean"].values for r, c in col.items()}
    lo = {r: t[f"{c}_CI_Low"].values for r, c in col.items()}
    hi = {r: t[f"{c}_CI_High"].values for r, c in col.items()}
    return names, means, lo, hi


def baseline() -> None:
    names, means, lo, hi = _baseline_ci()
    pf.plot_baseline_panels(means, names, ["ECM", "GF", "TF", "CYT"], OUT / "BSL_ECM_GF.png", lo, hi)
    pf.plot_baseline_panels(means, names, ["MET", "ION", "OXI", "CSF"], OUT / "BSL_MET_OXI.png", lo, hi)


def falsification() -> None:
    df = pd.read_csv(TABLES / "falsification_rules.csv")
    names, means, lo, hi = _baseline_ci()
    idx = [names.index(n) for n in df["node"]]
    pf.plot_falsification_profiles(
        df["node"], df["class"], df["delta"], df["ci_lower"], df["ci_upper"],
        means["Normal"][idx], lo["Normal"][idx], hi["Normal"][idx],
        means["Hyper"][idx], lo["Hyper"][idx], hi["Hyper"][idx],
        OUT / "NP_MT_FALS_RECT.png",
    )
    pf.plot_falsification_forest(df["node"], df["class"], df["delta"], df["ci_lower"],
                                 df["ci_upper"], df["outcome"], OUT / "NP_MT_FALS.png")


def rescue() -> None:
    rankings = pd.read_csv(TABLES / "true_rescue_rankings.csv")
    for tag, cat in pf.MODULE_CATEGORY.items():
        top = rankings[rankings["category"] == cat].sort_values("rank").head(20)
        pf.plot_true_rescue(list(top["rescue_strategy"]), top["rescue_percent"].values,
                            pf.RESCUE_TITLES[tag], OUT / f"Top20_{tag}_TrueRescue.png")

    z = np.load(REPS / "rescue_screen.npz")
    names = list(z["node_names"])
    labels = list(z["labels"])
    base = z["hyper_states"]                       # (R, N)
    pert = z["perturbed_states"]                   # (C, R, N)
    delta = pert - base[None]
    delta_mean = delta.mean(axis=1).T              # (N, C)
    delta_sd = delta.std(axis=1, ddof=1).T

    # Node-resolved responses of the two highest-ranked interventions per module.
    for fname, spec in pf.NODE_RESOLVED_FIGURES.items():
        panels = []
        for tag, title in spec:
            cat = pf.MODULE_CATEGORY[tag]
            nodes = [n for n in RESCUE_GROUPS[cat] if n in names]
            rows = [names.index(n) for n in nodes]
            top2 = rankings[rankings["category"] == cat].sort_values("rank").head(2)
            series = []
            for lab, pct in zip(top2["rescue_strategy"], top2["rescue_percent"]):
                c = labels.index(lab)
                series.append({"label": lab, "percent": float(pct),
                               "mean": delta_mean[rows, c], "sd": delta_sd[rows, c]})
            panels.append({"title": title, "nodes": nodes, "series": series})
        pf.plot_node_resolved_pairs(panels, OUT / f"{fname}.png")

    # Δ heatmaps with permutation-test BH-FDR stars (heatmap_with_sig).
    rng = np.random.RandomState(PERMUTATION_SEED)
    pmat = np.column_stack([
        permutation_pvalues(pert[c], base, n_perm=N_PERM, rng=rng)
        for c in range(pert.shape[0])
    ])
    qmat = bh_fdr(pmat.ravel()).reshape(pmat.shape)
    heat = OUT / "heatmaps"
    for tag, cat in pf.MODULE_CATEGORY.items():
        nodes = [n for n in RESCUE_GROUPS[cat] if n in names]
        rows = [names.index(n) for n in nodes]
        title = pf.RESCUE_TITLES[tag]
        for kind, cols in (("all", slice(None)), ("singles", slice(0, N_SINGLES)),
                           ("duals", slice(N_SINGLES, None))):
            sub = {"all": title, "singles": f"{title}: Single-node interventions",
                   "duals": f"{title}: Dual-node interventions"}[kind]
            pf.plot_rescue_heatmap(delta_mean[rows][:, cols], qmat[rows][:, cols],
                                   labels[cols], nodes, sub,
                                   heat / f"heatmap_{kind}_{tag}.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for step in (topology, baseline, falsification, rescue):
        print(f"[paper figures] {step.__name__}")
        step()
    print(f"[paper figures] written to {OUT}")


if __name__ == "__main__":
    main()
