"""Draw every computational figure of the manuscript, named as in the paper.

Reads the outputs of run_topology / run_baseline / run_transitions /
run_falsification / run_rescue (no re-simulation) and writes
results/figures/paper/<ManuscriptName>.png.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm import paper_figures as pf
from np_mt_rnm.network import load_network
from np_mt_rnm.statistics import bh_fdr, permutation_pvalues

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
REPS = RESULTS / "replicates"
OUT = RESULTS / "figures" / "paper"

# rng(PERMUTATION_SEED, 'twister') in RESCUE_NEW4_1_final.m, BASE_SEED + 4000.
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


def baseline() -> None:
    means = {}
    names = None
    for regime in ("Hypo", "Normal", "Hyper"):
        z = np.load(REPS / f"baseline_{regime.lower()}.npz")
        means[regime] = z["steady_states"].mean(axis=0)
        names = list(z["node_names"])
    pf.plot_baseline_panels(means, names, ["ECM", "GF", "TF", "CYT"], OUT / "BSL_ECM_GF.png")
    pf.plot_baseline_panels(means, names, ["MET", "ION", "OXI", "CSF"], OUT / "BSL_MET_OXI.png")


def transitions() -> None:
    h2n = pd.read_csv(TABLES / "transition_hypo_to_normal.csv").drop(columns="step")
    n2h = pd.read_csv(TABLES / "transition_normal_to_hyper.csv").drop(columns="step")
    names = list(h2n.columns)
    for group in ("ECM", "TF", "GF", "CYT", "OXI", "CSF"):
        pf.plot_transition_figure(h2n.values, n2h.values, names, group,
                                  OUT / f"{group}_transition.png")
    pf.plot_representative_paths(h2n.values, n2h.values, names,
                                 OUT / "Representative_transition_paths.png")


def falsification() -> None:
    df = pd.read_csv(TABLES / "falsification_rules.csv")
    pf.plot_falsification(df["node"], df["class"], df["delta"], df["ci_lower"],
                          df["ci_upper"], OUT / "NP_MT_FALS.png")


def rescue() -> None:
    df = pd.read_csv(TABLES / "true_rescue_rankings.csv")
    for group, cat in pf.MODULE_CATEGORY.items():
        top = df[df["category"] == cat].sort_values("rank")
        pf.plot_true_rescue(
            [pf.matlab_label(s) for s in top["rescue_strategy"]],
            top["rescue_percent"].values,
            pf.RESCUE_TITLES[group],
            OUT / f"{group}_rescue.png",
        )

    z = np.load(REPS / "rescue_screen.npz")
    names = list(z["node_names"])
    labels = [pf.matlab_label(s) for s in z["labels"]]
    base = z["hyper_states"]                       # (R, N)
    pert = z["perturbed_states"]                   # (C, R, N)
    delta = pert - base[None]
    delta_mean = delta.mean(axis=1).T              # (N, C)
    delta_sd = delta.std(axis=1, ddof=1).T
    rng = np.random.RandomState(PERMUTATION_SEED)
    pmat = np.column_stack([
        permutation_pvalues(pert[c], base, n_perm=N_PERM, rng=rng)
        for c in range(pert.shape[0])
    ])
    qmat = bh_fdr(pmat.ravel()).reshape(pmat.shape)
    for group in pf.MODULE_CATEGORY:
        pf.plot_node_resolved(group, names, labels, delta_mean, delta_sd, qmat,
                              OUT / f"{group}_rescue1.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for step in (topology, baseline, transitions, falsification, rescue):
        print(f"[paper figures] {step.__name__}")
        step()
    print(f"[paper figures] written to {OUT}")


if __name__ == "__main__":
    main()
