"""Reproduce paper Figs 4 and 5 (baseline steady-state activations).

Outputs:
  results/figures/fig4_baseline_ecm_gf_tf_cyto.png
  results/figures/fig5_baseline_metab_ion_ox_apop.png
  results/tables/baseline_summary.csv
  results/replicates/baseline_{hypo,normal,hyper}.npz
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.figures import plot_baseline_by_categories
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "MT_PRIMARY4_1.xlsx"
RESULTS = ROOT / "results"

N_REPS = 100
SEED = 20260420


def main() -> None:
    net = load_network(DATA)
    ensembles = {}
    for regime, preset in REGIME_PRESETS.items():
        print(f"[baseline] running {regime} ({N_REPS} replicates)...")
        ens = run_replicates(net, regime=preset, n_reps=N_REPS, seed=SEED, n_jobs=-1)
        assert ens.all_converged, f"{regime}: {(~ens.converged).sum()} non-converged replicates"
        np.savez_compressed(
            RESULTS / "replicates" / f"baseline_{regime.lower()}.npz",
            steady_states=ens.steady_states,
            node_names=np.array(ens.node_names),
            regime=np.array([f"{k}={v}" for k, v in preset.items()]),
        )
        ensembles[regime] = ens

    # Figure 4 — ECM + growth factors + transcription factors + cytokines/chemokines/proteases
    plot_baseline_by_categories(
        ensembles,
        categories=["ecm_matrix", "growth_factor", "transcription_factor", "cytokines_chemokines_proteases"],
        out_path=RESULTS / "figures" / "fig4_baseline_ecm_gf_tf_cyto.png",
        figsize=(14, 10),
    )
    # Figure 5 — metabolic + ion/Ca + oxidative + cell-fate
    plot_baseline_by_categories(
        ensembles,
        categories=["metabolic", "ion_channel", "oxidative_proteostasis", "cell_fate"],
        out_path=RESULTS / "figures" / "fig5_baseline_metab_ion_ox_apop.png",
        figsize=(14, 10),
    )

    # Summary CSV: one row per (regime, node) with mean + std.
    rows = []
    for regime, ens in ensembles.items():
        for i, name in enumerate(ens.node_names):
            rows.append({
                "regime": regime,
                "node": name,
                "mean_activation": float(ens.mean()[i]),
                "std_activation": float(ens.std()[i]),
            })
    pd.DataFrame(rows).to_csv(RESULTS / "tables" / "baseline_summary.csv", index=False)

    print("[baseline] done")


if __name__ == "__main__":
    main()
