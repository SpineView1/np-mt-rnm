"""Baseline steady states under Hypo / Normal / Hyper loading (manuscript Figs. 3-4).

Ports the baseline part of legacy/NP_MT_FALSIFY4_1_39.m:
  - 100 replicates per regime from random initial states (MATLAB-twister seeds);
  - 95% percentile bootstrap CIs of each node's mean activation
    (bootstrap_mean_CI, 10,000 resamples per condition);
  - Kruskal-Wallis across the three regimes and pairwise Wilcoxon rank-sum
    tests, BH-FDR adjusted across nodes within each test family
    (baseline_nonparametric_stats). The fixed loading inputs are not tested.

Outputs:
  results/replicates/baseline_{hypo,normal,hyper}.npz
  results/tables/baseline_summary.csv             (mean, SD per regime and node)
  results/tables/baseline_means_95ci.csv          (Baseline_Means_95CI.csv)
  results/tables/baseline_nonparametric_statistics.csv
  results/tables/baseline_statistical_summary.txt
The manuscript-styled panels are drawn by scripts/run_paper_figures.py.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import MECHANICAL_INPUTS, REGIME_PRESETS, REGIME_SEEDS, run_replicates
from np_mt_rnm.statistics import baseline_nonparametric_stats, bootstrap_mean_ci

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "MT_PRIMARY4_1.xlsx"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"

N_REPS = 100
BASE_NBOOT = 10_000
BASE_ALPHA = 0.05
BASE_STATS_ALPHA = 0.05
CI_SEED = 20261127


def main() -> None:
    net = load_network(DATA)
    ensembles = {}
    for regime, preset in REGIME_PRESETS.items():
        print(f"[baseline] running {regime} ({N_REPS} replicates)...")
        ens = run_replicates(
            net, regime=preset, n_reps=N_REPS, seed=REGIME_SEEDS[regime], n_jobs=-1
        )
        assert ens.all_converged, f"{regime}: {(~ens.converged).sum()} non-converged replicates"
        np.savez_compressed(
            RESULTS / "replicates" / f"baseline_{regime.lower()}.npz",
            steady_states=ens.steady_states,
            node_names=np.array(ens.node_names),
            regime=np.array([f"{k}={v}" for k, v in preset.items()]),
        )
        ensembles[regime] = ens

    names = list(net.node_names)

    rows = []
    for regime, ens in ensembles.items():
        for i, name in enumerate(ens.node_names):
            rows.append({
                "regime": regime,
                "node": name,
                "mean_activation": float(ens.mean()[i]),
                "std_activation": float(ens.std()[i]),
            })
    pd.DataFrame(rows).to_csv(TABLES / "baseline_summary.csv", index=False)

    # 95% bootstrap CIs of the mean (one generator, regimes in MATLAB order).
    rng = np.random.default_rng(CI_SEED)
    ci = {}
    for regime in ("Hypo", "Normal", "Hyper"):
        ci[regime] = bootstrap_mean_ci(ensembles[regime].steady_states, BASE_NBOOT,
                                       BASE_ALPHA, rng=rng)
    col = {"Hypo": "Hypo", "Normal": "NL", "Hyper": "HL"}
    table = {"Node": names}
    for regime in ("Hypo", "Normal", "Hyper"):
        table[f"{col[regime]}_Mean"] = ensembles[regime].mean()
        table[f"{col[regime]}_SD"] = ensembles[regime].std()
        table[f"{col[regime]}_CI_Low"] = ci[regime][0]
        table[f"{col[regime]}_CI_High"] = ci[regime][1]
    pd.DataFrame(table).to_csv(TABLES / "baseline_means_95ci.csv", index=False)

    stats = baseline_nonparametric_stats(
        ensembles["Hypo"].steady_states, ensembles["Normal"].steady_states,
        ensembles["Hyper"].steady_states, names, exclude=MECHANICAL_INPUTS,
    )
    stats.to_csv(TABLES / "baseline_nonparametric_statistics.csv", index=False)
    n_kw = int((stats.q_KW_3group < BASE_STATS_ALPHA).sum())
    n_nh = int((stats.q_NL_HL < BASE_STATS_ALPHA).sum())
    n_elig = int(np.isfinite(stats.p_KW_3group).sum())
    summary = "\n".join([
        "Baseline comparisons of RNM outcomes across randomized initial conditions",
        f"Mean activation error bars: percentile bootstrap {100 * (1 - BASE_ALPHA):.1f}% CIs, "
        f"{BASE_NBOOT} resamples/condition",
        "Kruskal-Wallis: three conditions (Hypo, NL, HL) per node.",
        "Pairwise Wilcoxon rank-sum: Hypo-NL, Hypo-HL, NL-HL per node.",
        "BH-FDR correction across nodes within each comparison family.",
        f"Network nodes: {len(stats)}; eligible tested nodes: {n_elig} (fixed loading inputs excluded).",
        f"KW: q<{BASE_STATS_ALPHA:.3f} for {n_kw} nodes.",
        f"NL-HL: q<{BASE_STATS_ALPHA:.3f} for {n_nh} nodes.",
        "These p/q values describe model-output variation among initializations, "
        "NOT biological replicates.",
        "",
    ])
    (TABLES / "baseline_statistical_summary.txt").write_text(summary)
    print(summary)
    print("[baseline] done")


if __name__ == "__main__":
    main()
