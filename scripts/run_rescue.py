"""Reproduce paper Figs 10, 11, 12 (rescue screen + Supp S7–S9)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.categories import (
    CATEGORY_LABELS,
    CATEGORY_ORDER,
    NODE_CATEGORIES,
    nodes_in_category,
)
from np_mt_rnm.figures import plot_rescue_category, plot_true_rescue_ranking
from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import (
    DEFAULT_TOP_N,
    enumerate_perturbations,
    mean_abs_displacement,
    rank_group_rescue,
    run_perturbation,
)
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates

ROOT = Path(__file__).resolve().parents[1]
N_REPS = 100


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    perts = list(enumerate_perturbations())

    # MATLAB solves SS_hyper_all and SS_norm_all once, then reuses them for
    # every perturbation and for the TRUE-rescue distances. Do the same.
    print("[rescue] Hyper baseline (shared by all 35 perturbations) ...")
    hyper_baseline = run_replicates(
        net, REGIME_PRESETS["Hyper"], n_reps=N_REPS, seed=REGIME_SEEDS["Hyper"], n_jobs=-1
    )
    print("[rescue] Normal baseline (TRUE-rescue reference) ...")
    normal_baseline = run_replicates(
        net, REGIME_PRESETS["Normal"], n_reps=N_REPS, seed=REGIME_SEEDS["Normal"], n_jobs=-1
    )

    results = []
    for i, p in enumerate(perts):
        print(f"[rescue] {i+1}/{len(perts)}: {p.label}")
        result = run_perturbation(
            net,
            anabolic_up=p.anabolic_up,
            catabolic_down=p.catabolic_down,
            n_reps=N_REPS,
            n_jobs=-1,
            baseline=hyper_baseline,
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

    # ---- TRUE Hyper -> Normal rescue rankings, per biological group --------
    # Ports analyze_group_true_rescue in legacy/RESCUE_NEW4_1_final.m.
    hyper_mean = hyper_baseline.mean()
    normal_mean = normal_baseline.mean()
    # FINAL: (n_nodes, n_strategies), mean perturbed state per strategy.
    final_states = np.column_stack([r.perturbed_mean for r in results])
    labels = [r.perturbation.label for r in results]

    tables = ROOT / "results" / "tables"
    figs_rescue = figs / "true_rescue"
    figs_rescue.mkdir(parents=True, exist_ok=True)

    ranking_rows = []
    for cat in CATEGORY_ORDER:
        if cat == "other":
            continue
        ranking = rank_group_rescue(
            group=CATEGORY_LABELS[cat],
            group_nodes=nodes_in_category(cat),
            node_names=results[0].node_names,
            hyper_mean=hyper_mean,
            normal_mean=normal_mean,
            final_states=final_states,
            strategy_labels=labels,
            top_n=DEFAULT_TOP_N,
        )
        if ranking.missing_nodes:
            print(f"[rescue]   {cat}: nodes not in network: {ranking.missing_nodes}")

        per_group = pd.DataFrame({
            "rank": range(1, len(ranking.strategy_labels) + 1),
            "rescue_strategy": ranking.strategy_labels,
            "rescue_percent": ranking.rescue_percent,
            "distance_to_normal": ranking.distance_to_normal,
        })
        per_group.to_csv(tables / f"true_rescue_top{DEFAULT_TOP_N}_{cat}.csv", index=False)

        plot_true_rescue_ranking(
            ranking, figs_rescue / f"true_rescue_top{DEFAULT_TOP_N}_{cat}.png"
        )

        for rank, (lab, pct, dist) in enumerate(
            zip(ranking.strategy_labels, ranking.rescue_percent,
                ranking.distance_to_normal), start=1
        ):
            ranking_rows.append({
                "category": cat,
                "group": CATEGORY_LABELS[cat],
                "rank": rank,
                "rescue_strategy": lab,
                "rescue_percent": pct,
                "distance_to_normal": dist,
                "distance_hyper_to_normal": ranking.distance_hyper_to_normal,
            })
        best = ranking.strategy_labels[0]
        print(f"[rescue]   {cat}: best = {best} ({ranking.rescue_percent[0]:.1f}%)")

    pd.DataFrame(ranking_rows).to_csv(tables / "true_rescue_rankings.csv", index=False)

    # Replicate-level screen for the node-resolved manuscript figures
    # (*_rescue1.png): paired deltas need every replicate, not just means.
    np.savez_compressed(
        ROOT / "results" / "replicates" / "rescue_screen.npz",
        node_names=np.array(results[0].node_names),
        labels=np.array(labels),
        hyper_states=hyper_baseline.steady_states,
        normal_states=normal_baseline.steady_states,
        perturbed_states=np.stack([r.perturbed_states for r in results]),
    )

    print("[rescue] done")


if __name__ == "__main__":
    main()
