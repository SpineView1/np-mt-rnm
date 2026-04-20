"""Matplotlib figure code tuned to match the paper's visual style.

Colors, fonts, and panel arrangement reproduce the published figures
(Workineh & Noailly 2026). Each plotting function reads from a
ReplicateEnsemble (or equivalent data) and writes a PNG to a given path.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from np_mt_rnm.categories import CATEGORY_LABELS, nodes_in_category

if TYPE_CHECKING:
    from np_mt_rnm.simulation import ReplicateEnsemble

# Paper's regime color palette (eye-match to Figs 4–5).
REGIME_COLORS = {
    "Hypo":   "#3f72af",   # blue
    "Normal": "#1f9d55",   # green
    "Hyper":  "#c74343",   # red
}

FIG_DPI = 300


def _apply_paper_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def plot_baseline_by_categories(
    ensembles: dict[str, "ReplicateEnsemble"],
    categories: list[str],
    out_path: Path,
    panel_titles: dict[str, str] | None = None,
    figsize: tuple[float, float] = (14, 10),
    n_cols: int = 2,
) -> None:
    """Grouped bar plot: one panel per category, bars per node, grouped by regime.

    Reproduces paper Figs 4 (A–D) and 5 (A–D) when called with the
    appropriate subset of categories.
    """
    _apply_paper_style()
    n_panels = len(categories)
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, constrained_layout=True)
    axes = np.array(axes).reshape(-1)

    any_ens = next(iter(ensembles.values()))
    node_names = any_ens.node_names

    for ax, category in zip(axes, categories):
        # Use only nodes that are in this category AND exist in the network.
        cat_nodes = [n for n in nodes_in_category(category) if n in node_names]
        if not cat_nodes:
            ax.set_visible(False)
            continue
        x = np.arange(len(cat_nodes))
        width = 0.28
        regimes = list(ensembles.keys())
        for k, regime in enumerate(regimes):
            ens = ensembles[regime]
            means = np.array([ens.mean()[node_names.index(n)] for n in cat_nodes])
            stds = np.array([ens.std()[node_names.index(n)] for n in cat_nodes])
            ax.bar(
                x + (k - (len(regimes) - 1) / 2) * width,
                means,
                width,
                yerr=stds,
                label=regime,
                color=REGIME_COLORS[regime],
                edgecolor="black",
                linewidth=0.5,
                error_kw=dict(lw=0.6, capsize=2),
            )
        ax.set_xticks(x)
        ax.set_xticklabels(cat_nodes, rotation=60, ha="right")
        title = (panel_titles or {}).get(category, CATEGORY_LABELS[category])
        ax.set_title(title)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Activation")
    for ax in axes[n_panels:]:
        ax.set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.99, 1.0), frameon=False)
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_fig3_topology(
    metrics,   # TopologyMetrics — forward ref to avoid import cycle
    out_path: Path,
    top_n: int = 20,
) -> None:
    """Three-panel horizontal bar plot: signed out-degree, betweenness, harmonic closeness."""
    _apply_paper_style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 8), constrained_layout=True)

    def _barh(ax, scores: dict[str, float], title: str, xlabel: str, color: str):
        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
        names = [k for k, _ in ranked][::-1]   # reverse so highest is at top
        vals = [v for _, v in ranked][::-1]
        ax.barh(names, vals, color=color, edgecolor="black", linewidth=0.5)
        ax.set_title(title)
        ax.set_xlabel(xlabel)

    _barh(axes[0], metrics.signed_out_degree, "Signed out-degree", "edges", "#555")
    _barh(axes[1], metrics.betweenness, "Betweenness centrality", "", "#7a3b9f")
    _barh(axes[2], metrics.harmonic_closeness, "Harmonic closeness", "", "#2a7f62")
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
