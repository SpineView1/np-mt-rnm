"""Quick-look figures for analyses outside the manuscript's figure set.

The manuscript's own figures are drawn by np_mt_rnm.paper_figures. This module
keeps the transition-path heatmaps and trajectories (scripts/run_transitions.py)
and the regime colour palette used by the webapp bundle.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from np_mt_rnm.categories import nodes_in_category

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


def plot_transition_heatmap(
    result,   # TransitionResult (forward ref)
    category: str,
    out_path: Path,
    panel_title: str,
) -> None:
    """Single-panel heatmap: rows = nodes in category, cols = path steps."""
    _apply_paper_style()
    nodes = [n for n in nodes_in_category(category) if n in result.node_names]
    if not nodes:
        raise ValueError(f"no nodes in category {category!r} found in result")
    idx = [result.node_names.index(n) for n in nodes]
    data = result.means[:, idx].T   # (n_nodes, n_steps)

    fig, ax = plt.subplots(
        figsize=(6, max(3.0, 0.3 * len(nodes) + 2.0)), constrained_layout=True
    )
    im = ax.imshow(data, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(nodes)))
    ax.set_yticklabels(nodes)
    ax.set_xticks(np.arange(len(result.path)))
    ax.set_xticklabels([f"S{k+1}" for k in range(len(result.path))])
    ax.set_xlabel("Path step")
    ax.set_title(panel_title)
    fig.colorbar(im, ax=ax, label="Mean activation")
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_transition_trajectories(
    hypo_to_normal,   # TransitionResult
    normal_to_hyper,  # TransitionResult
    nodes: list[str],
    out_path: Path,
    anabolic: set[str] | None = None,
) -> None:
    """Two-panel line plot (Fig 8 A/B) of representative node trajectories.

    Anabolic nodes are drawn with solid lines; all others dashed. This
    matches the paper's convention (Fig 8 caption).
    """
    _apply_paper_style()
    anabolic = anabolic or set()
    fig, axes = plt.subplots(
        1, 2, figsize=(13, 5), constrained_layout=True, sharey=True
    )

    for ax, result, title in [
        (axes[0], hypo_to_normal, "Hypo → Normal"),
        (axes[1], normal_to_hyper, "Normal → Hyper"),
    ]:
        steps = np.arange(len(result.path))
        for node in nodes:
            if node not in result.node_names:
                continue
            j = result.node_names.index(node)
            style = "-" if node in anabolic else "--"
            ax.plot(steps, result.means[:, j], style, label=node, lw=1.6)
        ax.set_xticks(steps)
        ax.set_xticklabels([f"S{k+1}" for k in steps])
        ax.set_xlabel("Path step")
        ax.set_title(title)
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel("Mean activation")
    axes[1].legend(
        loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=7, frameon=False
    )
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
