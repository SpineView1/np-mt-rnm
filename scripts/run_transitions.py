"""Constrained loading-path transitions (Hypo→Normal, Normal→Hyper).

Ports legacy/NP_MT_RNM_FSA4_1.m. This analysis appeared in earlier manuscript
versions; the current manuscript does not report it. It is kept as an
additional analysis: per-step tables in results/tables/transition_*.csv and
heatmaps in results/figures/transitions/.
"""
from __future__ import annotations

from pathlib import Path

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
# NP_MT_RNM_FSA4_1.m draws its initial states with an unseeded rand inside
# parfor, so its exact ensembles cannot be regenerated even in MATLAB. We seed
# every replicate the way RESCUE_NEW4_1_48_targets_HQ.m does (MATLAB twister, one seed
# per replicate); step k of a path uses seeds base + 1000*k + r.
SEED_HYPO_TO_NORMAL = 10001
SEED_NORMAL_TO_HYPER = 20001


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    print("[transitions] running Hypo→Normal path ...")
    h2n = run_transition_path(
        net, HYPO_TO_NORMAL_PATH, n_reps=N_REPS, seed=SEED_HYPO_TO_NORMAL, n_jobs=-1
    )
    print("[transitions] running Normal→Hyper path ...")
    n2h = run_transition_path(
        net, NORMAL_TO_HYPER_PATH, n_reps=N_REPS, seed=SEED_NORMAL_TO_HYPER, n_jobs=-1
    )

    figs = ROOT / "results" / "figures" / "transitions"
    figs.mkdir(parents=True, exist_ok=True)
    for cat in ("ecm_matrix", "transcription_factor", "growth_factor",
                "cytokines_chemokines_proteases", "oxidative_proteostasis", "cell_fate"):
        plot_transition_heatmap(h2n, cat, figs / f"transition_{cat}_hypo_to_normal.png",
                                "Hypo → Normal")
        plot_transition_heatmap(n2h, cat, figs / f"transition_{cat}_normal_to_hyper.png",
                                "Normal → Hyper")

    # Representative trajectories (anabolic solid / catabolic dashed)
    nodes = [
        "COL2A1", "ACAN", "TIMP3",
        "COL1A1", "COL10A1",
        "SOX9", "NF-κB", "ROS",
        "IGF1", "VEGF",
        "Bcl2", "CASP3",
    ]
    anabolic = {"COL2A1", "ACAN", "TIMP3", "SOX9", "IGF1", "Bcl2"}
    plot_transition_trajectories(
        h2n, n2h, nodes, figs / "transition_trajectories.png", anabolic=anabolic
    )

    # Dump per-step mean activation as CSV (for downloads tab + webapp consumption).
    for name, result in [("hypo_to_normal", h2n), ("normal_to_hyper", n2h)]:
        df = pd.DataFrame(result.means, columns=result.node_names)
        df.insert(0, "step", range(len(result.path)))
        df.to_csv(ROOT / "results" / "tables" / f"transition_{name}.csv", index=False)

    print("[transitions] done")


if __name__ == "__main__":
    main()
