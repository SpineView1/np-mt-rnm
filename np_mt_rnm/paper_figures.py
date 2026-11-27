"""Manuscript figures, drawn to match the published PNGs one for one.

Each function writes one file named exactly like the manuscript figure
(BSL_ECM_GF.png, NP_MT_FALS.png, Top20_ECM_TrueRescue.png, ...). Styling
mirrors the MATLAB plotting helpers that produced them:

  - grouped baseline bars with 95% bootstrap CIs:
        plot_grouped_bars_with_CI                    (NP_MT_FALSIFY4_1_39.m)
  - falsification NL/HL profiles and sign grid:
        plot_falsification_AB_v12                    (NP_MT_FALSIFY4_1_39.m)
  - falsification forest plot (vertical CIs, zero reference only):
        "Forest plot: Normal vs Hyper"               (NP_MT_FALSIFY4_1_39.m)
  - Top-20 rescue rankings:
        analyze_group_true_rescue                    (RESCUE_NEW4_1_48_targets_HQ.m)
  - node-resolved responses of the two top-ranked interventions per module
        (manuscript Supplementary Figures; paired Δ, mean ± SD)
  - rescue Δ heatmaps with FDR stars (all 48 / singles / duals):
        heatmap_with_sig                             (RESCUE_NEW4_1_48_targets_HQ.m)
  - transition heatmaps (not part of the current manuscript):
        draw_custom_heatmap                          (NP_MT_RNM_FSA4_1.m)
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

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
# Baseline panels use NP_MT_FALSIFY4_1_39.m group_categories (the metabolic
# group ends in MitD and the ion-channel group includes AQP1/AQP5, as in
# BSL_MET_OXI.png). Rescue figures use the rescue script's own readout groups
# (np_mt_rnm.rescue.RESCUE_GROUPS).
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

# trueRescueTitles in RESCUE_NEW4_1_48_targets_HQ.m, for the eight modules the
# manuscript reports, keyed by the manuscript's figure tag (Top20_<tag>_TrueRescue).
RESCUE_TITLES = {
    "ECM": "ECM Anabolism & Phenotype Markers",
    "GF": "Growth Factors",
    "TF": "Transcription Factors",
    "CYT": "Cytokines, Chemokines, Proteases & Others",
    "MET": "Metabolic & Related",
    "ION": "Ion Channels & Related",
    "OXI": "Oxidative-Stress Defense & Proteostasis",
    "CS": "Cell Survival, Apoptosis & Mitophagy/DNA-Damage",
}

# Port category id for each manuscript rescue module.
MODULE_CATEGORY = {
    "ECM": "ecm_matrix",
    "GF": "growth_factor",
    "TF": "transcription_factor",
    "CYT": "cytokines_chemokines_proteases",
    "MET": "metabolic",
    "ION": "ion_channel",
    "OXI": "oxidative_proteostasis",
    "CS": "cell_fate",
}

# Panel titles of the node-resolved supplementary figures, in panel order.
NODE_RESOLVED_FIGURES = {
    "NodeResolved_ECM_GF_TF_CYT": [
        ("ECM", "ECM anabolism and phenotype markers"),
        ("GF", "Growth factors"),
        ("TF", "Transcription factors"),
        ("CYT", "Cytokines, chemokines and proteases"),
    ],
    "NodeResolved_MET_ION_SURV_OXI": [
        ("MET", "Metabolic and related"),
        ("ION", "Ion channels and related"),
        ("CS", "Cell survival, apoptosis and mitophagy"),
        ("OXI", "Oxidative-stress defense and proteostasis"),
    ],
}

ANAB_UP = ("SOX9", "PPARγ", "HIF-1α", "NRF2", "SIRT1", "FOXO")
CATAB_KD = ("RhoA-E", "PIEZO1", "ROS", "P38", "mTORC1", "NF-κB")

# plot_falsification_AB_v12 colours.
FALS_GREEN = (0.00, 0.69, 0.14)
FALS_RED = (0.92, 0.075, 0.095)
FALS_DARKGREEN = (0.00, 0.34, 0.06)
FALS_DARKRED = (0.58, 0.00, 0.04)
FALS_GRAY = (0.78, 0.80, 0.82)
FALS_DARK = (0.14, 0.14, 0.14)
FALS_EDGE = (0.45, 0.45, 0.45)


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
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
def _grouped_bars(ax, labels, Y, title, letter, lo=None, hi=None):
    n = len(labels)
    x = np.arange(1, n + 1)
    w = 0.86 / 3
    for k, col in enumerate((COL_HYPO, COL_NORM, COL_HYPER)):
        xk = x + (k - 1) * w
        ax.bar(xk, Y[:, k], w * 0.92, color=col, edgecolor="none")
        if lo is not None and hi is not None:
            # Asymmetric 95% percentile-bootstrap CIs of the mean.
            err = np.vstack([np.maximum(0, Y[:, k] - lo[:, k]), np.maximum(0, hi[:, k] - Y[:, k])])
            ax.errorbar(xk, Y[:, k], yerr=err, fmt="none", ecolor=(0.1, 0.1, 0.1),
                        elinewidth=0.75, capsize=2.5)
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
                         groups: Sequence[str], out_path: Path,
                         ci_low: dict[str, np.ndarray] | None = None,
                         ci_high: dict[str, np.ndarray] | None = None) -> None:
    """2×2 grouped Hypo/Normal/Hyper bars; dict keys Hypo/Normal/Hyper.

    With `ci_low`/`ci_high`, each bar carries its 95% bootstrap CI of the mean,
    as plot_grouped_bars_with_CI in NP_MT_FALSIFY4_1_39.m draws.
    """
    _style()
    idx = {n: i for i, n in enumerate(node_names)}
    regimes = ("Hypo", "Normal", "Hyper")
    fig, axes = plt.subplots(2, 2, figsize=(16, 8.6),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.08})
    for ax, key, letter in zip(axes.flat, groups, "ABCD"):
        title, nodes = GROUPS[key]
        nodes = [n for n in nodes if n in idx]
        rows = [idx[n] for n in nodes]
        Y = np.column_stack([np.asarray(means[r])[rows] for r in regimes])
        lo = hi = None
        if ci_low is not None and ci_high is not None:
            lo = np.column_stack([np.asarray(ci_low[r])[rows] for r in regimes])
            hi = np.column_stack([np.asarray(ci_high[r])[rows] for r in regimes])
        _grouped_bars(ax, nodes, Y, title, letter, lo, hi)
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
    fig, axes = plt.subplots(1, 2, figsize=(20, height), gridspec_kw={"wspace": 0.12})
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
                             gridspec_kw={"wspace": 0.0})
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
# NP_MT_FALS_RECT.png / NP_MT_FALS.png
# =============================================================================
def plot_falsification_profiles(nodes: Sequence[str], classes: Sequence[str],
                                delta: Sequence[float], ci_lo: Sequence[float],
                                ci_hi: Sequence[float],
                                mean_nl: Sequence[float], nl_lo: Sequence[float],
                                nl_hi: Sequence[float],
                                mean_hl: Sequence[float], hl_lo: Sequence[float],
                                hl_hi: Sequence[float], out_path: Path) -> None:
    """NP_MT_FALS_RECT.png: (A) parallel NL/HL mean activations with 95% CIs;
    (B) sign grid, upper-left triangle = RNM prediction, lower-right = experiment.

    Ports plot_falsification_AB_v12 (Figure_A_NL_HL_Parallel_v12 and
    Figure_B_Vertical_Signs_v12) of NP_MT_FALSIFY4_1_39.m: benchmarks sorted by
    Δ descending; green = higher under that loading condition, red = lower,
    grey = simulated CI includes zero.
    """
    _style()
    nodes = np.array(nodes)
    classes = np.array(classes)
    delta = np.asarray(delta, dtype=float)
    order = np.argsort(-delta, kind="stable")
    arr = lambda v: np.asarray(v, dtype=float)[order]
    nodes, classes = nodes[order], classes[order]
    lo, hi = arr(ci_lo), arr(ci_hi)
    nl, nll, nlh = arr(mean_nl), arr(nl_lo), arr(nl_hi)
    hl, hll, hlh = arr(mean_hl), arr(hl_lo), arr(hl_hi)
    n = len(nodes)

    fig = plt.figure(figsize=(20, 0.24 * n + 2.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.2], wspace=0.12)

    # --- A: parallel bars -------------------------------------------------
    ax = fig.add_subplot(gs[0])
    x_hl = 1.24
    y = np.arange(1, n + 1)
    for i in range(1, n + 1, 2):
        ax.axhspan(i - 0.5, i + 0.5, color=(0.979, 0.979, 0.979), zorder=0, lw=0)
    for off, vals, l, h, col, dcol in [(0.0, nl, nll, nlh, FALS_GREEN, FALS_DARKGREEN),
                                       (x_hl, hl, hll, hlh, FALS_RED, FALS_DARKRED)]:
        for t in np.arange(0, 1.01, 0.2):
            ax.axvline(off + t, color=(0.90, 0.90, 0.90), lw=0.42, zorder=0)
        ax.barh(y, vals, 0.66, left=off, color=col, edgecolor="none", zorder=2)
        ax.errorbar(off + vals, y, xerr=np.vstack([np.maximum(0, vals - l), np.maximum(0, h - vals)]),
                    fmt="none", ecolor=dcol, elinewidth=1.1, capsize=2.5, zorder=3)
        ax.text(off + 0.5, 0.15, {0.0: "NL State"}.get(off, "HL State"), ha="center",
                va="bottom", fontsize=20, fontweight="bold", color=FALS_DARK)
    ax.set_ylim(n + 0.5, -0.9)
    ax.set_xlim(-0.02, x_hl + 1.04)
    ticks = list(np.arange(0, 1.01, 0.2)) + list(x_hl + np.arange(0, 1.01, 0.2))
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"] * 2, fontsize=13)
    ax.set_yticks(y)
    ax.set_yticklabels(nodes, fontsize=10.5)
    ax.tick_params(direction="out", length=2)
    ax.set_xlabel("Mean steady-state activation (0–1)", fontsize=15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.06, 1.0, "(A)", transform=ax.transAxes, fontsize=22, fontweight="bold",
            ha="right", va="bottom")

    # --- B: sign grid, split into two blocks as in the manuscript ---------
    ax = fig.add_subplot(gs[1])
    exp_dir = np.where(classes == "anabolic", 1, -1)
    sim_dir = np.zeros(n, dtype=int)
    sim_dir[lo > 0] = 1
    sim_dir[hi < 0] = -1

    def colour(d, code):
        # code 2 = NL column, code 1 = HL column.
        if d == 0:
            return FALS_GRAY
        return FALS_GREEN if (d == 1 and code == 2) or (d == -1 and code == 1) else FALS_RED

    half = int(np.ceil(n / 2))
    blocks = [(0, range(half)), (3.4, range(half, n))]
    col_w, gap = 0.95, 0.18
    for ox, rows in blocks:
        bounds = [(ox + 0.56, ox + 0.56 + col_w), (ox + 0.56 + col_w + gap, ox + 0.56 + 2 * col_w + gap)]
        for r, i in enumerate(rows):
            yt, yb = r + 0.5, r + 1.5
            for col, (xl, xr) in enumerate(bounds):
                code = 2 - col
                ax.fill([xl, xr, xl], [yt, yt, yb], color=colour(sim_dir[i], code), lw=0)
                ax.fill([xr, xr, xl], [yt, yb, yb], color=colour(exp_dir[i], code), lw=0)
                ax.plot([xr, xl], [yt, yb], color=FALS_EDGE, lw=0.76)
            ax.text(ox + 0.48, r + 1, nodes[i], ha="right", va="center", fontsize=9.6)
        m = len(rows)
        for xl, xr in bounds:
            ax.plot([xl, xl], [0.5, m + 0.5], color=FALS_EDGE, lw=1.65)
            ax.plot([xr, xr], [0.5, m + 0.5], color=FALS_EDGE, lw=1.65)
            for j in range(m + 1):
                ax.plot([xl, xr], [j + 0.5, j + 0.5], color=FALS_EDGE, lw=0.9)
        for (xl, xr), lab in zip(bounds, ("NL State", "HL State")):
            ax.text((xl + xr) / 2, 0.2, lab, ha="center", va="bottom", fontsize=13)
    support = int(np.sum(sim_dir == exp_dir))
    opposite = int(np.sum(sim_dir == -exp_dir))
    uncertain = int(np.sum(sim_dir == 0))
    ky = half + 1.6
    for k, (c, lab) in enumerate([(FALS_GREEN, "Higher under loading"),
                                  (FALS_RED, "Lower under loading"),
                                  (FALS_GRAY, "Simulated CI includes zero")]):
        x0 = 0.56 + 2.25 * k
        ax.fill([x0, x0 + 0.3, x0 + 0.3, x0], [ky - 0.35, ky - 0.35, ky + 0.35, ky + 0.35],
                color=c, ec=(0.5, 0.5, 0.5), lw=0.5)
        ax.text(x0 + 0.4, ky, lab, va="center", fontsize=8.5, fontweight="normal")
    ax.text(0.56, ky + 1.2, "Upper-left triangle: RNM predicted; lower-right: experiment-derived",
            fontsize=8.5, va="center", fontweight="normal")
    ax.text(0.56, ky + 2.3, f"Supported: {support}   Contradicted: {opposite}   "
            f"Inconclusive: {uncertain}", fontsize=10, va="center")
    ax.set_ylim(ky + 3.0, -0.4)
    ax.set_xlim(-0.6, 6.9)
    ax.axis("off")
    ax.text(0.0, 1.0, "(B)", transform=ax.transAxes, fontsize=22, fontweight="bold",
            ha="right", va="bottom")
    _save(fig, out_path)


def plot_falsification_forest(nodes: Sequence[str], classes: Sequence[str],
                              delta: Sequence[float], ci_lo: Sequence[float],
                              ci_hi: Sequence[float], outcomes: Sequence[str],
                              out_path: Path) -> None:
    """NP_MT_FALS.png: tested nodes on x, Δ with vertical 95% bootstrap CIs.

    Green / red markers: upward ("+") / downward ("-") literature response;
    squares mark contradictory predictions. Zero is the only reference line
    (no tolerance band), as in NP_MT_FALSIFY4_1_39.m.
    """
    _style()
    nodes = np.array(nodes)
    classes = np.array(classes)
    outcomes = np.array(outcomes)
    delta = np.asarray(delta, dtype=float)
    lo, hi = np.asarray(ci_lo, dtype=float), np.asarray(ci_hi, dtype=float)
    order = np.argsort(-delta, kind="stable")
    fig, ax = plt.subplots(figsize=(20, 8.8))
    x = np.arange(1, len(order) + 1)
    grey = (0.55, 0.55, 0.55)
    for xk, k in zip(x, order):
        ax.plot([xk, xk], [lo[k], hi[k]], color=grey, lw=1.5)
        for v in (lo[k], hi[k]):
            ax.plot([xk - 0.16, xk + 0.16], [v, v], color=grey, lw=1.2)
    for xk, k in zip(x, order):
        col = COL_ANAB if classes[k] == "anabolic" else COL_CATA
        mk = "s" if outcomes[k] == "Contradictory" else "o"
        ax.plot(xk, delta[k], mk, mfc=col, mec="k", ms=11 if mk == "o" else 10, mew=1.2,
                zorder=3)
    ax.axhline(0, color="k", lw=1.1)
    ax.set_xticks(x)
    ax.set_xticklabels(nodes[order], rotation=90, fontsize=14)
    ax.set_xlim(0.4, len(order) + 0.6)
    ax.set_ylim(-1, 1.02)
    ax.set_yticks(np.arange(-1, 1.01, 0.2))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.10g}" if abs(v) > 1e-9 else "0"))
    ax.tick_params(labelsize=14)
    ax.set_xlabel("Tested nodes", fontsize=17)
    ax.set_ylabel("Mean difference (Normal - Hyper)", fontsize=17)
    ax.set_title("Falsification forest plot (Normal - Hyper, 95% CI)", fontsize=17)
    ax.grid(True, ls="--", alpha=0.20, color="k")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(1.5)
    _save(fig, out_path)


# =============================================================================
# Top20_<tag>_TrueRescue.png (analyze_group_true_rescue)
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


def plot_node_resolved_pairs(panels: Sequence[dict], out_path: Path) -> None:
    """NodeResolved_*.png: four modules, the two top-ranked interventions each.

    Each entry of `panels` has keys ``title``, ``nodes`` (module nodes, plotted
    top to bottom), and ``series``: two dicts with ``label`` (MATLAB short
    label), ``percent`` (module rescue %), ``mean`` and ``sd`` (paired Δ =
    intervention − HL per node, across the 100 matched replicates).
    Blue circles = first-ranked, orange squares = second-ranked; horizontal
    error bars = SD.
    """
    _style()
    styles = [{"marker": "o", "color": (0.13, 0.36, 0.69), "mfc": (0.13, 0.36, 0.69), "ms": 9},
              {"marker": "s", "color": (0.85, 0.37, 0.12), "mfc": "white", "ms": 7}]
    heights = [max(4.5, 0.42 * len(p["nodes"]) + 2.0) for p in panels]
    row_h = [max(heights[0], heights[1]), max(heights[2], heights[3])]
    fig = plt.figure(figsize=(20, sum(row_h) + 1.5))
    gs = fig.add_gridspec(2, 2, height_ratios=row_h, hspace=0.42, wspace=0.18)
    for k, (p, letter) in enumerate(zip(panels, "ABCD")):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        nodes = list(p["nodes"])
        y = np.arange(len(nodes))
        for s_i, (ser, st) in enumerate(zip(p["series"], styles)):
            off = -0.16 if s_i == 0 else 0.16
            m = np.asarray(ser["mean"], dtype=float)
            sd = np.asarray(ser["sd"], dtype=float)
            ax.errorbar(m, y + off, xerr=sd, fmt=st["marker"], color=st["color"],
                        mfc=st["mfc"], mec=st["color"], ms=st["ms"], mew=1.6,
                        elinewidth=1.6, capsize=4,
                        label=f"{s_i + 1}  {ser['label']}  ({ser['percent']:.1f}%)")
        ax.axvline(0, color=(0.45, 0.47, 0.5), lw=1.6, zorder=0)
        ax.set_yticks(y)
        ax.set_yticklabels(nodes, fontsize=13, fontweight="normal")
        ax.set_ylim(len(nodes) - 0.4, -0.6)
        ax.tick_params(axis="x", labelsize=13)
        ax.xaxis.grid(True, color=(0.85, 0.85, 0.9), lw=1.2)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Paired change from HL (activation, a.u.)", fontsize=15)
        ax.set_title(f"({letter}) {p['title']}", fontsize=17)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False,
                  fontsize=13, ncol=1)
    _save(fig, out_path)


def plot_rescue_heatmap(delta: np.ndarray, qmat: np.ndarray, clamp_labels: Sequence[str],
                        nodes: Sequence[str], title: str, out_path: Path,
                        alpha_q: float = 0.05) -> None:
    """heatmap_with_sig of RESCUE_NEW4_1_48_targets_HQ.m.

    `delta`, `qmat`: (n_nodes, n_interventions). Signed blue-white-red map with
    symmetric limits; stars = BH-FDR q < 0.05 / 0.01 / 0.001 (white on strongly
    coloured tiles).
    """
    _style()
    neg = np.column_stack([np.linspace(0.12, 1, 128), np.linspace(0.35, 1, 128),
                           np.linspace(0.75, 1, 128)])
    pos = np.column_stack([np.linspace(1, 0.78, 128), np.linspace(1, 0.18, 128),
                           np.linspace(1, 0.18, 128)])
    cmap = matplotlib.colors.ListedColormap(np.vstack([neg, pos]))
    n_int, n_nodes = delta.shape[1], delta.shape[0]
    fig, ax = plt.subplots(figsize=(max(10, 2.5 + 0.32 * n_int), max(5, 2.8 + 0.36 * n_nodes)))
    vmax = max(np.nanmax(np.abs(delta)), 1e-9)
    im = ax.imshow(delta, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto",
                   interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("Change in activation (perturbed - HL)", fontsize=11)
    for (i, j), q in np.ndenumerate(qmat):
        t = _stars(q, alpha_q)
        if t:
            ax.text(j, i, t, ha="center", va="center", fontsize=7, fontweight="bold",
                    color="w" if abs(delta[i, j]) > 0.65 * vmax else "k")
    ax.set_xticks(range(n_int))
    ax.set_xticklabels(clamp_labels, rotation=55, ha="right", rotation_mode="anchor",
                       fontsize=9)
    ax.set_yticks(range(n_nodes))
    ax.set_yticklabels(nodes, fontsize=10)
    ax.set_title(title, fontsize=13)
    _save(fig, out_path)
