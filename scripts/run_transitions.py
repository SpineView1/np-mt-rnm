"""Reproduce paper Figs 6, 7, 8 (+ Supp S3–S6)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.figures import (
    plot_transition_heatmap,
    plot_transition_trajectories,
)
from np_mt_rnm.network import load_network
from np_mt_rnm.transitions import (
    HYPO_TO_NORMAL_PATH,
    NORMAL_TO_HYPER_PATH,
    run_transition_path,
)

ROOT = Path(__file__).resolve().parents[1]
N_REPS = 100
SEED = 20260420


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    print("[transitions] running Hypo→Normal path ...")
    h2n = run_transition_path(
        net, HYPO_TO_NORMAL_PATH, n_reps=N_REPS, seed=SEED, n_jobs=-1
    )
    print("[transitions] running Normal→Hyper path ...")
    n2h = run_transition_path(
        net, NORMAL_TO_HYPER_PATH, n_reps=N_REPS, seed=SEED + 10000, n_jobs=-1
    )

    figs = ROOT / "results" / "figures"
    figs_supp = figs / "supp"
    figs_supp.mkdir(parents=True, exist_ok=True)

    # Figure 6 — ECM heatmaps (A: Hypo→Normal, B: Normal→Hyper)
    plot_transition_heatmap(
        h2n, "ecm_matrix", figs / "fig6_transition_ecm_A.png", "Hypo → Normal"
    )
    plot_transition_heatmap(
        n2h, "ecm_matrix", figs / "fig6_transition_ecm_B.png", "Normal → Hyper"
    )

    # Figure 7 — Transcription factor heatmaps
    plot_transition_heatmap(
        h2n,
        "transcription_factor",
        figs / "fig7_transition_tf_A.png",
        "Hypo → Normal",
    )
    plot_transition_heatmap(
        n2h,
        "transcription_factor",
        figs / "fig7_transition_tf_B.png",
        "Normal → Hyper",
    )

    # Supplementary S3–S6 — remaining major categories
    for cat, sfig in [
        ("growth_factor", "S3"),
        ("cytokines_chemokines_proteases", "S4"),
        ("oxidative_proteostasis", "S5"),
        ("cell_fate", "S6"),
    ]:
        plot_transition_heatmap(
            h2n, cat, figs_supp / f"{sfig}_h2n.png", "Hypo → Normal"
        )
        plot_transition_heatmap(
            n2h, cat, figs_supp / f"{sfig}_n2h.png", "Normal → Hyper"
        )

    # Figure 8 — representative trajectories (anabolic solid / catabolic dashed)
    nodes = [
        "COL2A1", "ACAN", "TIMP3",
        "COL1A1", "COL10A1",
        "SOX9", "NF-κB", "ROS",
        "IGF1", "VEGF",
        "Bcl2", "CASP3",
    ]
    anabolic = {"COL2A1", "ACAN", "TIMP3", "SOX9", "IGF1", "Bcl2"}
    plot_transition_trajectories(
        h2n, n2h, nodes, figs / "fig8_transition_trajectories.png", anabolic=anabolic
    )

    # Dump per-step mean activation as CSV (for downloads tab + webapp consumption).
    for name, result in [("hypo_to_normal", h2n), ("normal_to_hyper", n2h)]:
        df = pd.DataFrame(result.means, columns=result.node_names)
        df.insert(0, "step", range(len(result.path)))
        df.to_csv(ROOT / "results" / "tables" / f"transition_{name}.csv", index=False)

    print("[transitions] done")


if __name__ == "__main__":
    main()
