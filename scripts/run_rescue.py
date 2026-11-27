"""Hyper→Normal rescue screen: 12 single + 36 dual node clamps (manuscript Figs. 6-7).

Ports legacy/RESCUE_NEW4_1_48_targets_HQ.m. Writes:
  results/tables/rescue_intervention_design_48.csv  (Rescue_Intervention_Design_48.csv)
  results/tables/rescue_per_node_deltas.csv         (mean paired Δ, 48 × nodes)
  results/tables/rescue_per_node_sd.csv             (SD of paired Δ)
  results/tables/rescue_summary.csv                 (mean |Δ| per module)
  results/tables/true_rescue_rankings.csv           (all 48 ranked, every module)
  results/tables/true_rescue_top20_<module>.csv
  results/tables/rescue_module_top.csv              (best intervention per module)
  results/tables/rescue_phenotype_panels.csv        (anabolic / catabolic panel means)
  results/replicates/rescue_screen.npz
Manuscript figures are drawn by scripts/run_paper_figures.py.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from np_mt_rnm.categories import CATEGORY_LABELS
from np_mt_rnm.network import load_network
from np_mt_rnm.rescue import (
    ANABOLIC_PANEL,
    CATABOLIC_PANEL,
    DEFAULT_TOP_N,
    N_PERTURBATIONS,
    N_SINGLES,
    PAPER_RESCUE_MODULES,
    RESCUE_GROUPS,
    enumerate_perturbations,
    mean_abs_displacement,
    rank_group_rescue,
    run_perturbation,
)
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
N_REPS = 100


def _panel_mean(names: list[str], vec: np.ndarray, panel) -> float:
    idx = [names.index(n) for n in panel if n in names]
    return float(np.mean(vec[idx]))


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    names = list(net.node_names)
    perts = list(enumerate_perturbations())
    assert len(perts) == N_PERTURBATIONS == 48

    design = pd.DataFrame({
        "InterventionID": range(1, len(perts) + 1),
        "InterventionType": [
            "Single suppression" if p.anabolic_up is None else
            "Single activation" if p.catabolic_down is None else
            "Dual activation + suppression" for p in perts
        ],
        "ShortLabel": [p.short_label for p in perts],
        "ClampAssignments": [
            ", ".join(([f"{p.anabolic_up}=1"] if p.anabolic_up else [])
                      + ([f"{p.catabolic_down}=0"] if p.catabolic_down else []))
            for p in perts
        ],
    })
    design.to_csv(TABLES / "rescue_intervention_design_48.csv", index=False)

    # MATLAB solves SS_hyper_all and SS_norm_all once, then reuses them for
    # every perturbation and for the TRUE-rescue distances. Do the same.
    print("[rescue] Hyper baseline (shared by all 48 perturbations) ...")
    hyper = run_replicates(
        net, REGIME_PRESETS["Hyper"], n_reps=N_REPS, seed=REGIME_SEEDS["Hyper"], n_jobs=-1
    )
    print("[rescue] Normal baseline (TRUE-rescue reference) ...")
    normal = run_replicates(
        net, REGIME_PRESETS["Normal"], n_reps=N_REPS, seed=REGIME_SEEDS["Normal"], n_jobs=-1
    )

    results = []
    for i, p in enumerate(perts):
        print(f"[rescue] {i + 1}/{len(perts)}: {p.short_label}")
        results.append(run_perturbation(
            net, anabolic_up=p.anabolic_up, catabolic_down=p.catabolic_down,
            n_reps=N_REPS, n_jobs=-1, baseline=hyper,
        ))
    labels = [r.perturbation.short_label for r in results]

    rows = []
    for r in results:
        row = {
            "perturbation": r.perturbation.short_label,
            "anabolic_up": r.perturbation.anabolic_up or "",
            "catabolic_down": r.perturbation.catabolic_down or "",
        }
        for cat in PAPER_RESCUE_MODULES:
            row[f"mean_abs_delta_{cat}"] = mean_abs_displacement(r, list(RESCUE_GROUPS[cat]))
        rows.append(row)
    pd.DataFrame(rows).to_csv(TABLES / "rescue_summary.csv", index=False)

    for fname, attr in (("rescue_per_node_deltas.csv", "mean_delta"),
                        ("rescue_per_node_sd.csv", "std_delta")):
        df = pd.DataFrame(np.stack([getattr(r, attr) for r in results]), columns=names)
        df.insert(0, "perturbation", labels)
        df.to_csv(TABLES / fname, index=False)

    # ---- TRUE Hyper -> Normal rescue, per readout group --------------------
    hyper_mean, normal_mean = hyper.mean(), normal.mean()
    final_states = np.column_stack([r.perturbed_mean for r in results])

    ranking_rows, top_rows = [], []
    for cat, group_nodes in RESCUE_GROUPS.items():
        ranking = rank_group_rescue(
            group=CATEGORY_LABELS[cat], group_nodes=list(group_nodes), node_names=names,
            hyper_mean=hyper_mean, normal_mean=normal_mean, final_states=final_states,
            strategy_labels=labels, top_n=len(labels),
        )
        if ranking.missing_nodes:
            print(f"[rescue]   {cat}: nodes not in network: {ranking.missing_nodes}")
        top = pd.DataFrame({
            "rank": range(1, DEFAULT_TOP_N + 1),
            "rescue_strategy": ranking.strategy_labels[:DEFAULT_TOP_N],
            "rescue_percent": ranking.rescue_percent[:DEFAULT_TOP_N],
            "distance_to_normal": ranking.distance_to_normal[:DEFAULT_TOP_N],
        })
        top.to_csv(TABLES / f"true_rescue_top{DEFAULT_TOP_N}_{cat}.csv", index=False)
        for rank, (lab, pct, dist) in enumerate(zip(
                ranking.strategy_labels, ranking.rescue_percent, ranking.distance_to_normal), 1):
            ranking_rows.append({
                "category": cat, "group": CATEGORY_LABELS[cat], "rank": rank,
                "rescue_strategy": lab, "rescue_percent": pct, "distance_to_normal": dist,
                "distance_hyper_to_normal": ranking.distance_hyper_to_normal,
                "in_manuscript": cat in PAPER_RESCUE_MODULES,
            })
        top_rows.append({
            "category": cat, "group": CATEGORY_LABELS[cat],
            "top_intervention": ranking.strategy_labels[0],
            "rescue_percent": ranking.rescue_percent[0],
            "in_manuscript": cat in PAPER_RESCUE_MODULES,
        })
        print(f"[rescue]   {cat}: best = {ranking.strategy_labels[0]} "
              f"({ranking.rescue_percent[0]:.1f}%)")
    pd.DataFrame(ranking_rows).to_csv(TABLES / "true_rescue_rankings.csv", index=False)
    (pd.DataFrame(top_rows).sort_values("rescue_percent", ascending=False)
       .to_csv(TABLES / "rescue_module_top.csv", index=False))

    # ---- Phenotype panels (Supplementary Table) ----------------------------
    pheno = [
        {"regime_or_intervention": "NL baseline", "type": "baseline",
         "anabolic_mean": _panel_mean(names, normal_mean, ANABOLIC_PANEL),
         "catabolic_mean": _panel_mean(names, normal_mean, CATABOLIC_PANEL)},
        {"regime_or_intervention": "HL baseline", "type": "baseline",
         "anabolic_mean": _panel_mean(names, hyper_mean, ANABOLIC_PANEL),
         "catabolic_mean": _panel_mean(names, hyper_mean, CATABOLIC_PANEL)},
    ]
    for j, lab in enumerate(labels):
        pheno.append({
            "regime_or_intervention": lab,
            "type": "single" if j < N_SINGLES else "dual",
            "anabolic_mean": _panel_mean(names, final_states[:, j], ANABOLIC_PANEL),
            "catabolic_mean": _panel_mean(names, final_states[:, j], CATABOLIC_PANEL),
        })
    pd.DataFrame(pheno).to_csv(TABLES / "rescue_phenotype_panels.csv", index=False)

    # Replicate-level screen for the node-resolved figures and FDR heatmaps.
    np.savez_compressed(
        ROOT / "results" / "replicates" / "rescue_screen.npz",
        node_names=np.array(names),
        labels=np.array(labels),
        hyper_states=hyper.steady_states,
        normal_states=normal.steady_states,
        perturbed_states=np.stack([r.perturbed_states for r in results]),
    )
    print("[rescue] done")


if __name__ == "__main__":
    main()
