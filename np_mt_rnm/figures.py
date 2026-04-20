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


def plot_fig9_forest(
    outcomes,   # list[FalsificationOutcome] — forward ref
    out_path: Path,
) -> None:
    """Two-panel Fig 9: (A) bar of Δ colored by class, (B) forest plot with 95% CIs."""
    _apply_paper_style()
    fig, axes = plt.subplots(
        1, 2, figsize=(14, max(6, 0.18 * len(outcomes))), constrained_layout=True
    )
    sorted_outs = sorted(outcomes, key=lambda o: o.delta)
    nodes = [o.rule.node for o in sorted_outs]
    deltas = [o.delta for o in sorted_outs]
    colors = [
        "#1f9d55" if o.rule.cls == "anabolic" else "#c74343" for o in sorted_outs
    ]

    axes[0].barh(nodes, deltas, color=colors, edgecolor="black", linewidth=0.3)
    axes[0].axvline(0, color="k", lw=0.6)
    axes[0].set_xlabel("Δ = Normal − Hyper")
    axes[0].set_title("A. Directional effects")

    for i, o in enumerate(sorted_outs):
        axes[1].plot([o.ci_lower, o.ci_upper], [i, i], color="k", lw=0.8)
        axes[1].plot(o.delta, i, "o", color=colors[i], markersize=4)
    axes[1].set_yticks(np.arange(len(sorted_outs)))
    axes[1].set_yticklabels(nodes)
    axes[1].axvline(0, color="k", lw=0.4, ls=":")
    axes[1].axvline(0.02, color="gray", lw=0.4, ls="--")
    axes[1].axvline(-0.02, color="gray", lw=0.4, ls="--")
    axes[1].set_xlabel("Δ with 95% bootstrap CI")
    axes[1].set_title("B. Bootstrap CIs")

    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_rescue_category(
    results,   # list[PerturbationResult] — forward ref
    category: str,
    out_path: Path,
    title: str,
) -> None:
    """4-panel figure reproducing Fig 10/11/12:
       A = heatmap of mean Δ across perturbations × nodes in category,
       B = ranked bar of mean|Δ| per perturbation,
       C = node-level Δ for the top-ranked combinatorial perturbation,
       D = node-level Δ for the strongest single-node perturbation.
    """
    from np_mt_rnm.rescue import mean_abs_displacement
    _apply_paper_style()

    node_names = results[0].node_names
    nodes = [n for n in nodes_in_category(category) if n in node_names]
    idx = [node_names.index(n) for n in nodes]
    if not idx:
        raise ValueError(f"no network nodes in category {category!r}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=True)

    # Panel A — heatmap of perturbations × nodes-in-category.
    M = np.array([r.mean_delta[idx] for r in results])
    labels = [r.perturbation.label for r in results]
    vmax = float(np.max(np.abs(M))) if np.max(np.abs(M)) > 0 else 1.0
    im = axes[0, 0].imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[0, 0].set_yticks(np.arange(len(labels)))
    axes[0, 0].set_yticklabels(labels, fontsize=7)
    axes[0, 0].set_xticks(np.arange(len(nodes)))
    axes[0, 0].set_xticklabels(nodes, rotation=60, ha="right", fontsize=7)
    axes[0, 0].set_title("A. Mean Δ (perturbed − Hyper baseline)")
    fig.colorbar(im, ax=axes[0, 0], label="Δ")

    # Panel B — ranked mean|Δ| per perturbation.
    scores = [mean_abs_displacement(r, nodes) for r in results]
    order = np.argsort(scores)[::-1]
    axes[0, 1].barh(
        [labels[i] for i in order][::-1],
        [scores[i] for i in order][::-1],
        color="#555",
        edgecolor="black",
        linewidth=0.4,
    )
    axes[0, 1].set_xlabel("Mean |Δ|")
    axes[0, 1].set_title(f"B. Ranked by mean|Δ| — {category}")

    # Panel C — highest-ranked combinatorial perturbation (anabolic↑ + catabolic↓).
    combos = [i for i in order if results[i].perturbation.anabolic_up and results[i].perturbation.catabolic_down]
    top_combo = results[combos[0]] if combos else results[order[0]]
    c_vals = [top_combo.mean_delta[node_names.index(n)] for n in nodes]
    c_errs = [top_combo.std_delta[node_names.index(n)] for n in nodes]
    axes[1, 0].bar(
        nodes,
        c_vals,
        yerr=c_errs,
        color=["#1f9d55" if v > 0 else "#c74343" for v in c_vals],
        edgecolor="black", linewidth=0.4, error_kw=dict(lw=0.6, capsize=2),
    )
    axes[1, 0].axhline(0, color="k", lw=0.5)
    axes[1, 0].set_xticklabels(nodes, rotation=60, ha="right")
    axes[1, 0].set_ylabel("Δ")
    axes[1, 0].set_title(f"C. Top combo: {top_combo.perturbation.label}")

    # Panel D — strongest single-node perturbation.
    singles = [i for i in order if not (results[i].perturbation.anabolic_up and results[i].perturbation.catabolic_down)]
    top_single = results[singles[0]] if singles else results[order[0]]
    s_vals = [top_single.mean_delta[node_names.index(n)] for n in nodes]
    s_errs = [top_single.std_delta[node_names.index(n)] for n in nodes]
    axes[1, 1].bar(
        nodes,
        s_vals,
        yerr=s_errs,
        color=["#1f9d55" if v > 0 else "#c74343" for v in s_vals],
        edgecolor="black", linewidth=0.4, error_kw=dict(lw=0.6, capsize=2),
    )
    axes[1, 1].axhline(0, color="k", lw=0.5)
    axes[1, 1].set_xticklabels(nodes, rotation=60, ha="right")
    axes[1, 1].set_ylabel("Δ")
    axes[1, 1].set_title(f"D. Top single: {top_single.perturbation.label}")

    fig.suptitle(title, fontsize=13)
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
