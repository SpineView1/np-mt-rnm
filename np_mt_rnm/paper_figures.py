"""Manuscript figures, drawn to match the published PNGs one for one.

Each function writes one file named exactly like the manuscript figure
(BSL_ECM_GF.png, ECM_rescue.png, ...). Styling mirrors the MATLAB plotting
helpers that produced them:

  - grouped baseline bars:    plot_grouped_bars_no_errors  (NP_MT_RNM_FALSIFY4_1.m)
  - falsification bar/forest: PART B figures               (NP_MT_RNM_FALSIFY4_1.m)
  - transition heatmaps:      draw_custom_heatmap          (NP_MT_RNM_FSA4_1.m)
  - Top-20 rescue rankings:   analyze_group_true_rescue    (RESCUE_NEW4_1_final.m)
  - node-resolved rescue:     heatmap_with_sig, bar_top,
                              bar_with_sd_and_sig          (RESCUE_NEW4_1_final.m)
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

DPI = 300

# --- MATLAB colours ---------------------------------------------------------
COL_HYPO = (0.0, 0.0, 1.0)
COL_NORM = (0.0, 0.60, 0.0)
COL_HYPER = (1.0, 0.0, 0.0)
COL_ANAB = (0.15, 0.65, 0.15)
COL_CATA = (0.85, 0.25, 0.25)
MATLAB_BLUE = (0.0, 0.447, 0.741)           # default bar() colour
MATLAB_ORANGE = (0.929, 0.694, 0.125)
MATLAB_GREEN = (0.466, 0.674, 0.188)
MATLAB_CYAN = (0.301, 0.745, 0.933)
ROYAL_BLUE = (0.255, 0.412, 0.882)

# MATLAB parula, 16 anchor points (R2014b+ default colormap).
_PARULA_ANCHORS = [
    (0.2422, 0.1504, 0.6603), (0.2810, 0.2178, 0.8291), (0.2795, 0.3074, 0.9418),
    (0.2365, 0.4057, 0.9883), (0.1707, 0.5003, 0.9600), (0.1468, 0.5862, 0.9025),
    (0.1155, 0.6623, 0.8605), (0.0280, 0.7178, 0.8075), (0.0711, 0.7586, 0.7232),
    (0.2209, 0.7929, 0.6158), (0.4161, 0.8132, 0.4722), (0.6377, 0.8031, 0.3233),
    (0.8148, 0.7716, 0.2142), (0.9608, 0.7471, 0.2203), (0.9789, 0.8506, 0.1726),
    (0.9769, 0.9839, 0.0805),
]
PARULA = LinearSegmentedColormap.from_list("parula", _PARULA_ANCHORS, N=256)

# --- Node groups as plotted in the manuscript ---------------------------------
# Baseline panels use NP_MT_RNM_FALSIFY4_1.m group_categories (the metabolic
# group ends in MitD and the ion-channel group includes AQP1/AQP5, as in
# BSL_MET_OXI.png). Transition and rescue figures use the same lists for the
# six modules they show.
GROUPS: dict[str, tuple[str, list[str]]] = {
    "ECM": ("ECM anabolism & phenotype markers",
            ["COL2A1", "COL1A1", "COL10A1", "ACAN", "TIMP3"]),
    "GF": ("Growth factors",
           ["TGFβ", "VEGF", "IGF1", "BMP2", "CCN2", "GDF5", "FGF2", "FGF18", "Wnt3a", "Wnt5a"]),
    "TF": ("Transcription Factors",
           ["CREB", "HIF-1α", "HIF-2α", "NF-κB", "AP-1", "FOXO", "SOX9", "NFAT", "RUNX2",
            "YAP/TAZ", "MRTF-A", "NRF2", "HSF1", "TonEBP", "ELK1", "PPARγ", "CITED2"]),
    "CYT": ("Cytokines, chemokines, proteases & others",
            ["TNF", "IL6", "IL1β", "IL8", "CCL2", "CXCL1", "CXCL3", "ADAMTS4/5",
             "MMP1", "MMP3", "MMP13"]),
    "MET": ("Metabolic & related",
            ["LKB1", "NAD+", "AMPK", "mTORC1", "mTORC2", "SIRT1", "PI3K-M", "PI3K-E",
             "PIP3-M", "PIP3-E", "PDK1-M", "PDK1-E", "AKT1-M", "AKT1-E", "GSK3B", "ULK1",
             "PTEN", "PLD2", "PGE2", "COX-2", "CAT", "GPX1", "SOD1", "SOD2", "HO-1",
             "PHD2", "VHL", "Rheb", "MitD"]),
    "ION": ("Ion channels & related",
            ["Ca2+os", "Ca2+su", "CaMKII", "PKC-E", "PKC-M", "PLCγ-M", "PLCγ-E", "CaN",
             "IP3", "PLA2", "AQP1", "AQP5"]),
    "OXI": ("Oxidative-stress defense & proteostasis",
            ["HO-1", "GPX1", "SOD1", "SOD2", "CAT", "HSP70", "HSP27", "ROS"]),
    "CSF": ("Cell survival, apoptosis & mitophagy/DNA-damage",
            ["Bcl2", "BAX", "CASP3", "CASP9", "BNIP3", "GADD45", "DRP1", "MOMP"]),
}

# trueRescueTitles in RESCUE_NEW4_1_final.m.
RESCUE_TITLES = {
    "ECM": "ECM Anabolism & Phenotype Markers",
    "GF": "Growth Factors",
    "TF": "Transcription Factors",
    "CYT": "Cytokines, Chemokines, Proteases & Others",
    "OXI": "Oxidative-Stress Defense & Proteostasis",
    "CSF": "Cell Survival, Apoptosis & Mitophagy/DNA-Damage",
}

# Port category id for each manuscript module (identical node lists).
MODULE_CATEGORY = {
    "ECM": "ecm_matrix",
    "GF": "growth_factor",
    "TF": "transcription_factor",
    "CYT": "cytokines_chemokines_proteases",
    "OXI": "oxidative_proteostasis",
    "CSF": "cell_fate",
}

ANAB_UP = ("SOX9", "PPARγ", "HIF-1α", "NRF2", "IκBα")

# Panel C of *_rescue1.png shows "the highest-ranked dual-node perturbation"
# by mean |Δ| (Supplementary S3). In four modules the manuscript figure follows
# that rule; for ECM and CYT it shows the second-ranked dual instead. These
# entries reproduce the manuscript as published.
MANUSCRIPT_PANEL_C = {"ECM": "SOX9-ROS", "CYT": "PPARγ-FAK-E"}
CATAB_KD = ("RhoA-E", "PIEZO1", "PI3K-E", "FAK-E", "ROS")


def _style() -> None:
    plt.rcParams.update({
        "font.family": "Arial",
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.unicode_minus": True,
        "mathtext.default": "regular",
    })


def _panel_letter(ax, letter: str, x: float = -0.02, y: float = 1.04, size: float = 22) -> None:
    ax.text(x, y, f"{letter})", transform=ax.transAxes, fontsize=size,
            fontweight="bold", ha="right", va="bottom")


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def matlab_label(port_label: str) -> str:
    """'SOX9↑ + FAK-E↓' -> 'SOX9-FAK-E' (ShortClampLabels in normalize_clamps)."""
    return port_label.replace("↑ + ", "-").replace("↑", "").replace("↓", "")


def _matlab_ticks(axis, lo: float, hi: float, nbins: int = 7, snap: bool = True):
    """MATLAB-like tick steps (1/2/5 x 10^k), 'g' labels, limits snapped to ticks."""
    loc = matplotlib.ticker.MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10])
    ticks = loc.tick_values(lo, hi)
    step = ticks[1] - ticks[0]
    if snap:
        lo, hi = np.floor(lo / step + 1e-9) * step, np.ceil(hi / step - 1e-9) * step
    first = np.ceil(lo / step - 1e-9) * step
    ticks = np.round(np.arange(first, hi + step * 1e-6, step), 10)
    axis.set_ticks(ticks)
    axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{(0.0 if abs(v) < 1e-12 else v):.10g}"))
    return lo, hi


def _fmt3g(v: float) -> str:
    """MATLAB sprintf('%.3g') (no leading zeros in the exponent padding)."""
    s = f"{v:.3g}"
    if "e" in s:
        mant, exp = s.split("e")
        sign = exp[0]
        digits = exp[1:].lstrip("0").rjust(2, "0")
        s = f"{mant}e{sign}{digits}"
    return s


# =============================================================================
# TOPO_STATS.png
# =============================================================================
def edge_list_node_order(net) -> list[str]:
    """Node order of a MATLAB digraph built from [activation edges; inhibition edges].

    Nodes enter in order of first appearance as an edge source (activations
    row by row, then inhibitions), then as a target. MATLAB's descending sort
    is stable, so this order breaks ties in TOPO_STATS panel A.
    """
    n = len(net.node_names)
    src, tgt = [], []
    for m in (net.mact, net.minh):
        for i in range(n):
            for j in np.nonzero(m[i])[0]:
                src.append(net.node_names[j])
                tgt.append(net.node_names[i])
    order = list(dict.fromkeys(src + tgt))
    return order + [x for x in net.node_names if x not in order]


def plot_topo_stats(node_names: Sequence[str], out_act: Sequence[int], out_inh: Sequence[int],
                    betweenness_raw: Sequence[float], harmonic_out_raw: Sequence[float],
                    out_path: Path, node_order: Sequence[str] | None = None,
                    exclude_from_closeness: Sequence[str] = ("Hypo", "NL", "HL"),
                    top_degree: int = 25, top_n: int = 20) -> None:
    """Signed out-degree split, raw betweenness, raw outgoing harmonic closeness.

    Arrays are reordered to `node_order` first so that ties sort as in MATLAB.
    The mechanical inputs are not ranked in panel C, as in the manuscript.
    """
    _style()
    names = np.array(node_names)
    if node_order is not None:
        perm = [list(names).index(x) for x in node_order]
        names = names[perm]
        out_act = np.asarray(out_act)[perm]
        out_inh = np.asarray(out_inh)[perm]
        betweenness_raw = np.asarray(betweenness_raw)[perm]
        harmonic_out_raw = np.asarray(harmonic_out_raw, dtype=float)[perm].copy()
    harmonic_out_raw = np.asarray(harmonic_out_raw, dtype=float).copy()
    harmonic_out_raw[np.isin(names, list(exclude_from_closeness))] = -np.inf
    act = np.asarray(out_act)
    inh = np.asarray(out_inh)
    tot = act + inh

    fig = plt.figure(figsize=(13.5, 10.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.42, wspace=0.16)

    ax = fig.add_subplot(gs[0, :])
    order = np.argsort(-tot, kind="stable")[:top_degree]
    x = np.arange(len(order))
    blue, red = (0.20, 0.51, 0.93), (0.90, 0.30, 0.30)
    ax.bar(x, act[order], 0.78, color=blue, edgecolor="k", linewidth=0.6, label="+")
    ax.bar(x, inh[order], 0.78, bottom=act[order], color=red, edgecolor="k", linewidth=0.6, label="−")
    for xi, k in zip(x, order):
        if act[k] > 0:
            ax.text(xi, act[k] / 2, str(act[k]), ha="center", va="center", color="white",
                    fontsize=20, fontweight="bold")
        if inh[k] > 0:
            ax.text(xi, act[k] + inh[k] / 2, str(inh[k]), ha="center", va="center",
                    color="white", fontsize=20, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(names[order], rotation=45, ha="right", rotation_mode="anchor", fontsize=20)
    ax.set_xlim(-0.8, len(order) - 0.2)
    ax.set_ylim(0, 5 * np.ceil((tot[order].max() + 0.5) / 5))
    ax.tick_params(axis="y", labelsize=20, width=2)
    ax.set_ylabel("Out-degree split", fontsize=28)
    ax.set_title("Top broadcasters (+/−)", fontsize=22)
    ax.grid(True, color=(0.8, 0.8, 0.8), linewidth=1.2)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_linewidth(2)
    h = [Patch(facecolor=red, edgecolor="k"), Patch(facecolor=blue, edgecolor="k")]
    leg = ax.legend(h, ["−", "+"], loc="upper right", fontsize=34, handlelength=2.6,
                    handleheight=1.5, frameon=True, fancybox=False, edgecolor="k")
    leg.get_frame().set_linewidth(3)
    _panel_letter(ax, "A", x=0.22, y=1.0, size=22)

    for col, (vals, title, ylab, letter) in enumerate([
        (np.asarray(betweenness_raw), "Top 20 by Betweenness", "Betweenness", "B"),
        (np.asarray(harmonic_out_raw), "Top 20 by HarmonicClosenessOut", "HarmonicClosenessOut", "C"),
    ]):
        ax = fig.add_subplot(gs[1, col])
        o = np.argsort(-vals, kind="stable")[:top_n]
        ax.bar(np.arange(len(o)), vals[o], 0.8, color=ROYAL_BLUE)
        ax.set_xticks(np.arange(len(o)))
        ax.set_xticklabels(names[o], rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
        ax.tick_params(axis="y", labelsize=9)
        ax.set_xlim(-0.8, len(o) - 0.2)
        ax.set_ylim(*_matlab_ticks(ax.yaxis, 0, vals[o].max(), nbins=10 if letter == "B" else 6))
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.grid(True, color=(0.85, 0.85, 0.85), linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        _panel_letter(ax, letter, x=0.27, y=1.035, size=22)

    _save(fig, out_path)


# =============================================================================
# BSL_ECM_GF.png / BSL_MET_OXI.png
# =============================================================================
def _grouped_bars(ax, labels, Y, title, letter):
    n = len(labels)
    x = np.arange(1, n + 1)
    w = 0.86 / 3
    for k, col in enumerate((COL_HYPO, COL_NORM, COL_HYPER)):
        ax.bar(x + (k - 1) * w, Y[:, k], w * 0.92, color=col, edgecolor="none")
    if n <= 8:
        rot = 0
    elif n <= 14:
        rot = 45
    else:
        rot = 60
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=rot, ha="center" if rot == 0 else "right",
                       rotation_mode="anchor", fontsize=12)
    ax.set_xlim(0.35, n + 0.65)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=12)
    ax.set_ylabel("Activation (a.u.)", fontsize=13)
    ax.tick_params(direction="out", width=1.2)
    ax.grid(True, color=(0.80, 0.80, 0.80), alpha=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.5)
    t = ax.set_title(title, fontsize=13, pad=8)
    ax.annotate(f"{letter})", xy=(0, 0.5), xycoords=t, xytext=(-8, 0),
                textcoords="offset points", ha="right", va="center",
                fontsize=19, fontweight="bold")


def plot_baseline_panels(means: dict[str, np.ndarray], node_names: Sequence[str],
                         groups: Sequence[str], out_path: Path) -> None:
    """2×2 grouped Hypo/Normal/Hyper bars; `means` keys Hypo/Normal/Hyper."""
    _style()
    idx = {n: i for i, n in enumerate(node_names)}
    fig, axes = plt.subplots(2, 2, figsize=(16, 8.6),
                             gridspec_kw=dict(hspace=0.55, wspace=0.08))
    for ax, key, letter in zip(axes.flat, groups, "ABCD"):
        title, nodes = GROUPS[key]
        nodes = [n for n in nodes if n in idx]
        Y = np.column_stack([np.asarray(means[r])[[idx[n] for n in nodes]]
                             for r in ("Hypo", "Normal", "Hyper")])
        _grouped_bars(ax, nodes, Y, title, letter)
    handles = [Patch(color=c) for c in (COL_HYPO, COL_NORM, COL_HYPER)]
    fig.legend(handles, ["Hypo", "Normal", "Hyper"], loc="upper center", ncol=3,
               frameon=False, fontsize=13, handlelength=4, bbox_to_anchor=(0.5, 0.995))
    _save(fig, out_path)


# =============================================================================
# *_transition.png
# =============================================================================
def _transition_heatmap(ax, Z, nodes, title, letter, fig):
    """draw_custom_heatmap: nodes × path steps, parula, YDir normal, %.3g cells."""
    nr, nc = Z.shape
    im = ax.imshow(Z, cmap=PARULA, aspect="auto", origin="lower",
                   extent=(0.5, nc + 0.5, 0.5, nr + 0.5), interpolation="nearest")
    vmin, vmax = np.nanmin(Z), np.nanmax(Z)
    if vmax <= vmin:
        vmax = vmin + 1e-12
    im.set_clim(vmin, vmax)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.012)
    _matlab_ticks(cb.ax.yaxis, vmin, vmax, nbins=10, snap=False)
    cb.ax.set_ylim(vmin, vmax)
    cb.ax.tick_params(labelsize=10)
    for t in cb.ax.get_yticklabels():
        t.set_fontweight("bold")
    for i in np.arange(0.5, nr + 1):
        ax.plot([0.5, nc + 0.5], [i, i], "k-", lw=0.5)
    for j in np.arange(0.5, nc + 1):
        ax.plot([j, j], [0.5, nr + 0.5], "k-", lw=0.5)
    midv = (vmin + vmax) / 2
    fs = 11 if nr <= 10 else 9
    for i in range(nr):
        for j in range(nc):
            v = Z[i, j]
            ax.text(j + 1, i + 1, _fmt3g(v), ha="center", va="center", fontsize=fs,
                    fontweight="bold", color="k" if v > midv else "w")
    ax.set_xticks(range(1, nc + 1))
    ax.set_yticks(range(1, nr + 1))
    ax.set_yticklabels(nodes)
    ax.tick_params(labelsize=11, width=1.2)
    ax.set_xlim(0.5, nc + 0.5)
    ax.set_ylim(0.5, nr + 0.5)
    ax.set_xlabel("Path step", fontsize=13)
    ax.set_title(title, fontsize=14)
    _panel_letter(ax, letter, x=0.08, y=1.0, size=20)


def plot_transition_figure(h2n: np.ndarray, n2h: np.ndarray, node_names: Sequence[str],
                           group: str, out_path: Path) -> None:
    """`h2n`/`n2h`: (n_steps, n_nodes) mean activations along each path."""
    _style()
    title, nodes = GROUPS[group]
    idx = {n: i for i, n in enumerate(node_names)}
    nodes = [n for n in nodes if n in idx]
    cols = [idx[n] for n in nodes]
    height = max(5.0, 0.62 * len(nodes) + 1.8)
    fig, axes = plt.subplots(1, 2, figsize=(20, height), gridspec_kw=dict(wspace=0.12))
    _transition_heatmap(axes[0], np.asarray(h2n)[:, cols].T, nodes,
                        f"{title} (Hypo to Normal)", "A", fig)
    _transition_heatmap(axes[1], np.asarray(n2h)[:, cols].T, nodes,
                        f"{title} (Normal to Hyper)", "B", fig)
    _save(fig, out_path)


# =============================================================================
# Representative_transition_paths.png
# =============================================================================
REPRESENTATIVE_NODES = [
    ("COL2A1", MATLAB_BLUE, "-"), ("ACAN", MATLAB_GREEN, "-"), ("TIMP3", MATLAB_ORANGE, "-"),
    ("COL1A1", MATLAB_BLUE, "--"), ("COL10A1", MATLAB_GREEN, "--"), ("SOX9", "k", "-"),
    ("NF-κB", MATLAB_ORANGE, "--"), ("ROS", "k", "--"), ("IGF1", "magenta", "-"),
    ("VEGF", "magenta", "--"), ("Bcl2", (0.0, 0.75, 1.0), "-"), ("CASP3", (0.0, 0.75, 1.0), "--"),
]


def plot_representative_paths(h2n: np.ndarray, n2h: np.ndarray, node_names: Sequence[str],
                              out_path: Path) -> None:
    _style()
    idx = {n: i for i, n in enumerate(node_names)}
    fig, axes = plt.subplots(1, 2, figsize=(20, 8.2), sharey=True,
                             gridspec_kw=dict(wspace=0.0))
    for ax, data, title, letter in [
        (axes[0], np.asarray(h2n), "Hypo to Normal transition path", "A"),
        (axes[1], np.asarray(n2h), "Normal to Hyper transition path", "B"),
    ]:
        steps = np.arange(1, data.shape[0] + 1)
        for node, col, ls in REPRESENTATIVE_NODES:
            ax.plot(steps, data[:, idx[node]], ls, color=col, lw=2.6, label=node)
        ax.set_xlim(1, data.shape[0])
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(1, data.shape[0] + 0.01, 0.5))
        ax.set_yticks(np.arange(0, 1.01, 0.1))
        ax.tick_params(labelsize=11)
        ax.set_xlabel("Path step", fontsize=13)
        ax.set_title(title, fontsize=14)
        ax.grid(True, color=(0.85, 0.85, 0.85))
        ax.text(0.22, 1.012, f"{letter})", transform=ax.transAxes, fontsize=22,
                fontweight="bold", ha="right", va="bottom")
    axes[0].set_ylabel("Steady-state activation", fontsize=13)
    axes[1].tick_params(axis="y", left=False)
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
    for ax in axes:
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)
    # Panel A's last tick sits on panel B's first; MATLAB shows only B's.
    axes[0].set_xticks(np.arange(1, h2n.shape[0] - 0.49, 0.5))
    handles, labels = axes[0].get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, 0.5),
                     fontsize=16, handlelength=3.0, frameon=True, fancybox=False,
                     edgecolor="k", framealpha=1.0)
    leg.get_frame().set_linewidth(1.5)
    pos = axes[0].get_position()
    for y0, y1 in ((pos.y1 - 0.02, pos.y1 + 0.03), (pos.y0 - 0.06, pos.y0 - 0.02)):
        fig.add_artist(plt.Line2D([pos.x1, pos.x1], [y0, y1], transform=fig.transFigure,
                                  color=(0.1, 0.35, 0.5), lw=4))
    _save(fig, out_path)


# =============================================================================
# NP_MT_FALS.png
# =============================================================================
def plot_falsification(nodes: Sequence[str], classes: Sequence[str], delta: Sequence[float],
                       ci_lo: Sequence[float], ci_hi: Sequence[float], out_path: Path,
                       ftol: float = 0.02) -> None:
    _style()
    nodes = np.array(nodes)
    classes = np.array(classes)
    delta = np.asarray(delta)
    lo, hi = np.asarray(ci_lo), np.asarray(ci_hi)
    order = np.argsort(-delta, kind="stable")
    cols = [COL_ANAB if c == "anabolic" else COL_CATA for c in classes[order]]

    fig = plt.figure(figsize=(20, 7.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.15, 1], wspace=0.13)

    ax = fig.add_subplot(gs[0])
    x = np.arange(1, len(order) + 1)
    ax.bar(x, delta[order], 0.92, color=cols, edgecolor="none")
    ax.axhline(0, color="k", lw=2.0)
    ax.axhline(ftol, color=(0.3, 0.9, 0.3), ls=":", lw=2.5)
    ax.axhline(-ftol, color=COL_CATA, ls=":", lw=2.5)
    ax.set_xticks(x)
    ax.set_xticklabels(nodes[order], rotation=60, ha="right", rotation_mode="anchor", fontsize=17)
    ax.set_xlim(0.3, len(order) + 0.7)
    ax.set_ylim(-1, 1)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticklabels(["-1", "-0.5", "0", "0.5", "1"], fontsize=18)
    ax.set_ylabel("Δ (Normal - Hyper)", fontsize=20)
    ax.set_title("Directional effects by node (Normal - Hyper, 95% CI)", fontsize=19)
    ax.grid(True, color=(0.80, 0.80, 0.80), alpha=0.5, lw=1.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(2)
    ax.text(0.08, 1.09, "A)", transform=ax.transAxes, fontsize=22, fontweight="bold")

    ax = fig.add_subplot(gs[1])
    y = np.arange(1, len(order) + 1)
    for yk, k in zip(y, order):
        ax.plot([lo[k], hi[k]], [yk, yk], color=(0.55, 0.55, 0.55), lw=2.2)
    for yk, k, c in zip(y, order, cols):
        ax.plot(delta[k], yk, "o", mfc=c, mec="k", ms=8, mew=0.9)
    ax.axvline(-ftol, color=COL_CATA, ls=":", lw=1.4)
    ax.axvline(ftol, color=(0.3, 0.9, 0.3), ls=":", lw=1.4)
    ax.set_yticks(y)
    ax.set_yticklabels(nodes[order], fontsize=10)
    ax.set_ylim(len(order) + 0.5, 0.5)
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xticklabels(["-1", "-0.5", "0", "0.5", "1"], fontsize=10)
    ax.set_xlabel("Mean difference (Normal - Hyper)", fontsize=11)
    ax.set_title("Falsification forest plot (Normal - Hyper, 95% CI)", fontsize=10.5)
    ax.grid(True, color=(0.80, 0.80, 0.80), alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.1, 1.07, "B)", transform=ax.transAxes, fontsize=22, fontweight="bold")
    _save(fig, out_path)


# =============================================================================
# *_rescue.png (analyze_group_true_rescue)
# =============================================================================
def plot_true_rescue(labels: Sequence[str], percents: Sequence[float], title: str,
                     out_path: Path, normal_ref_threshold: float = 80.0) -> None:
    """Top-N horizontal bars; strongest at the top. Labels in MATLAB format."""
    _style()
    best = np.asarray(percents, dtype=float)
    labels = list(labels)
    n_top = len(best)
    plot_score = best[::-1]
    plot_labels = labels[::-1]

    fig = plt.figure(figsize=(1900 / 96, 1250 / 96))
    ax = fig.add_axes([0.28, 0.11, 0.66, 0.80])
    colors = [(0.12, 0.40, 0.72) if s >= 0 else (0.78, 0.20, 0.20) for s in plot_score]
    ax.barh(np.arange(1, n_top + 1), plot_score, 0.68, color=colors, edgecolor="none")

    best_plot = plot_score.max()
    min_score = min(plot_score.min(), 0)
    upper_target = best_plot * 1.15 if best_plot > 0 else 10
    x_upper = max(np.ceil(upper_target / 10) * 10, 10)
    x_lower = np.floor((min_score * 1.15) / 10) * 10 if min_score < 0 else 0

    ax.axvline(0, color="k", ls="--", lw=2.5)
    if best_plot >= normal_ref_threshold:
        x_upper = max(x_upper, 105)
        ax.axvline(100, color="k", ls=":", lw=2.5)
        # xline label: MATLAB puts it right of the line, reading upward, at the top.
        ax.text(100 + 0.006 * (x_upper - x_lower), n_top + 0.6, "100% Normal restoration",
                rotation=90, ha="left", va="top", fontsize=20, fontweight="bold")
    ax.set_xlim(x_lower, x_upper)
    rng = max(x_upper - x_lower, 1)
    step = 10 if rng <= 60 else 20
    ax.set_xticks(np.arange(np.ceil(x_lower / step) * step,
                            np.floor(x_upper / step) * step + 0.01, step))
    ax.set_yticks(np.arange(1, n_top + 1))
    ax.set_yticklabels(plot_labels)
    ax.tick_params(labelsize=21, width=2, direction="out", length=6)
    ax.xaxis.grid(True, alpha=0.16, color="k", lw=1.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(2)
    ax.set_xlabel("Rescue toward Normal profile (%)", fontsize=26)
    ax.set_title(f"Top {n_top} Rescue Strategies — {title}\nHyper → Normal Profile Restoration",
                 fontsize=29)
    for k, s in enumerate(plot_score, start=1):
        xpos = s + 0.015 * rng if s >= 0 else s - 0.015 * rng
        ax.text(xpos, k, f"{s:.1f}%", fontsize=24, fontweight="bold",
                ha="left" if s >= 0 else "right", va="center", clip_on=False)
    ax.set_ylim(0.3, n_top + 0.7)
    _save(fig, out_path)


# =============================================================================
# *_rescue1.png (heatmap_with_sig + bar_top + bar_with_sd_and_sig)
# =============================================================================
def _stars(q: float, alpha: float = 0.05) -> str:
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < alpha:
        return "*"
    return ""


def _clamp_title(label: str) -> str:
    """format_clamp_title: '↑ SOX9, ↓ ROS' with -E/-M stripped."""
    # Split on the known node names: they contain '-' themselves (FAK-E, HIF-1α).
    parts = []
    rest = label
    for up in ANAB_UP:
        if rest == up or rest.startswith(up + "-"):
            parts.append("↑ " + up)
            rest = rest[len(up) + 1:]
            break
    for kd in CATAB_KD:
        if rest == kd:
            parts.append("↓ " + kd.removesuffix("-E").removesuffix("-M"))
    return ", ".join(parts)


def plot_node_resolved(group: str, node_names: Sequence[str], clamp_labels: Sequence[str],
                       delta_mean: np.ndarray, delta_sd: np.ndarray, qmat: np.ndarray,
                       out_path: Path, n_top: int = 10) -> None:
    """Four panels: A Δ heatmap with FDR stars, B top-10 by mean|Δ|, C top dual, D top single.

    `delta_mean`, `delta_sd`, `qmat`: (n_nodes_all, n_clamps); `clamp_labels`
    in MATLAB short form ('SOX9', 'SOX9-FAK-E') and enumeration order.
    """
    _style()
    title, nodes = GROUPS[group]
    idx = {n: i for i, n in enumerate(node_names)}
    nodes = [n for n in nodes if n in idx]
    rows = [idx[n] for n in nodes]
    D = delta_mean[rows]
    S = delta_sd[rows]
    Q = qmat[rows]
    labels = list(clamp_labels)

    fig = plt.figure(figsize=(20, 11.3))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.13, width_ratios=[1.06, 1])

    # A: heatmap_with_sig
    ax = fig.add_subplot(gs[0, 0])
    vmax = max(np.nanmax(np.abs(D)), 1e-9)
    im = ax.imshow(D, cmap=PARULA, vmin=-vmax, vmax=vmax, aspect="auto",
                   interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.008)
    _matlab_ticks(cb.ax.yaxis, -vmax, vmax, nbins=10, snap=False)
    cb.ax.set_ylim(-vmax, vmax)
    cb.ax.tick_params(labelsize=7)
    for (i, j), q in np.ndenumerate(Q):
        t = _stars(q)
        if t:
            ax.text(j, i, t, ha="center", va="center", fontsize=6, fontweight="bold")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=6.5)
    ax.set_yticks(range(len(nodes)))
    ax.set_yticklabels(nodes, fontsize=7)
    ax.tick_params(length=2)
    ax.set_title(title, fontsize=9)
    ax.text(0.13, 1.03, "A)", transform=ax.transAxes, fontsize=22, fontweight="bold")

    # B: bar_top
    scores = np.nanmean(np.abs(D), axis=0)
    top = np.argsort(-scores, kind="stable")[:n_top]
    ax = fig.add_subplot(gs[0, 1])
    ax.bar(np.arange(1, len(top) + 1), scores[top], 0.8, color=MATLAB_BLUE,
           edgecolor=(0.15, 0.15, 0.15), linewidth=0.5)
    ax.set_xticks(np.arange(1, len(top) + 1))
    ax.set_xticklabels([labels[k] for k in top], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=12)
    ax.set_xlim(0.2, len(top) + 0.8)
    ax.tick_params(labelsize=11)
    ax.set_ylim(*_matlab_ticks(ax.yaxis, 0, scores[top].max(), nbins=10))
    ax.set_ylabel("Mean |Δ|", fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.grid(True, color=(0.85, 0.85, 0.85))
    ax.set_axisbelow(True)
    ax.text(0.06, 1.03, "B)", transform=ax.transAxes, fontsize=22, fontweight="bold")

    # C/D: highest-ranked dual and single by the same mean|Δ| criterion.
    ranked = np.argsort(-scores, kind="stable")
    is_dual = np.array(["-" in lab and lab not in CATAB_KD and lab not in ANAB_UP
                        for lab in labels])
    c_dual = next(k for k in ranked if is_dual[k])
    if group in MANUSCRIPT_PANEL_C:
        c_dual = labels.index(MANUSCRIPT_PANEL_C[group])
    c_single = next(k for k in ranked if not is_dual[k])
    for spec, c, letter in [(gs[1, 0], c_dual, "C"), (gs[1, 1], c_single, "D")]:
        ax = fig.add_subplot(spec)
        m, s, q = D[:, c], S[:, c], Q[:, c]
        x = np.arange(1, len(nodes) + 1)
        ax.bar(x, m, 0.8, color=(0.25, 0.45, 0.85), edgecolor=(0.2, 0.2, 0.2), linewidth=0.4)
        ax.errorbar(x, m, yerr=s, fmt="none", ecolor="k", elinewidth=1.0, capsize=2)
        yoff = max(np.max(s), 1e-6) * 0.5 + 0.02
        for xi, mi, si, qi in zip(x, m, s, q):
            t = _stars(qi)
            if t:
                ax.text(xi, mi + si + yoff, t, ha="center", fontsize=8, fontweight="bold")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylim(*_matlab_ticks(ax.yaxis, min(0, np.min(m - s)), max(0, np.max(m + s)), nbins=6))
        ax.set_xticks(x)
        ax.set_xticklabels(nodes, rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
        ax.set_xlim(0, len(nodes) + 1)
        ax.tick_params(labelsize=9)
        ax.set_ylabel("Δ (mean ± SD)", fontsize=11)
        ax.set_title(f"{title}: {_clamp_title(labels[c])}", fontsize=12)
        ax.grid(True, color=(0.85, 0.85, 0.85))
        ax.set_axisbelow(True)
        ax.text(-0.04, 1.03, f"{letter})", transform=ax.transAxes, fontsize=22,
                fontweight="bold", ha="right")
    _save(fig, out_path)
