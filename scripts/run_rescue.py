"""Reproduce paper Figs 10, 11, 12 (rescue screen + Supp S7–S9)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.categories import NODE_CATEGORIES
from np_mt_rnm.figures import plot_rescue_category
from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import (
    enumerate_perturbations,
    mean_abs_displacement,
    run_perturbation,
)

ROOT = Path(__file__).resolve().parents[1]
N_REPS = 100
SEED = 20260420


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    perts = list(enumerate_perturbations())

    results = []
    for i, p in enumerate(perts):
        print(f"[rescue] {i+1}/{len(perts)}: {p.label}")
        result = run_perturbation(
            net,
            anabolic_up=p.anabolic_up,
            catabolic_down=p.catabolic_down,
            n_reps=N_REPS,
            seed=SEED + 1000 * i,
            n_jobs=-1,
        )
        results.append(result)

    figs = ROOT / "results" / "figures"
    figs_supp = figs / "supp"
    figs_supp.mkdir(parents=True, exist_ok=True)

    plot_rescue_category(results, "ecm_matrix", figs / "fig10_rescue_ecm.png",
                        "Figure 10. ECM-phenotype rescue")
    plot_rescue_category(results, "growth_factor", figs / "fig11_rescue_growth_factors.png",
                        "Figure 11. Growth-factor rescue")
    plot_rescue_category(results, "transcription_factor", figs / "fig12_rescue_transcription.png",
                        "Figure 12. Transcription-factor rescue")

    # Supplementary
    for cat, sfig in [
        ("cytokines_chemokines_proteases", "S7"),
        ("oxidative_proteostasis", "S8"),
        ("cell_fate", "S9"),
    ]:
        plot_rescue_category(
            results, cat, figs_supp / f"{sfig}_rescue.png",
            f"Supp {sfig}. Rescue — {cat}"
        )

    # Summary CSV — one row per perturbation with per-category mean|Δ| scores.
    def _nodes_for(cat):
        return [n for n, cs in NODE_CATEGORIES.items() if cat in cs]

    rows = []
    for r in results:
        rows.append({
            "perturbation": r.perturbation.label,
            "anabolic_up": r.perturbation.anabolic_up or "",
            "catabolic_down": r.perturbation.catabolic_down or "",
            "mean_abs_delta_ecm":
                mean_abs_displacement(r, _nodes_for("ecm_matrix")),
            "mean_abs_delta_tf":
                mean_abs_displacement(r, _nodes_for("transcription_factor")),
            "mean_abs_delta_gf":
                mean_abs_displacement(r, _nodes_for("growth_factor")),
            "mean_abs_delta_cytokines":
                mean_abs_displacement(r, _nodes_for("cytokines_chemokines_proteases")),
            "mean_abs_delta_oxidative":
                mean_abs_displacement(r, _nodes_for("oxidative_proteostasis")),
            "mean_abs_delta_cell_fate":
                mean_abs_displacement(r, _nodes_for("cell_fate")),
        })
    pd.DataFrame(rows).to_csv(ROOT / "results" / "tables" / "rescue_summary.csv", index=False)

    # Per-node deltas CSV (35 × 147).
    per_node = np.stack([r.mean_delta for r in results])
    df = pd.DataFrame(per_node, columns=results[0].node_names)
    df.insert(0, "perturbation", [r.perturbation.label for r in results])
    df.to_csv(ROOT / "results" / "tables" / "rescue_per_node_deltas.csv", index=False)

    print("[rescue] done")


if __name__ == "__main__":
    main()
