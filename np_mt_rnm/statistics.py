"""Statistical helpers used by falsification and rescue analyses.

Pure-numpy implementations for reproducibility and to avoid extra dependencies.
"""
from __future__ import annotations

import numpy as np


def bootstrap_diff_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int = 10_000,
    alpha: float = 0.05,
    rng: np.random.Generator | None = None,
) -> tuple[float, float, float]:
    """Bootstrap (alpha/2, 1-alpha/2) CI for x.mean() - y.mean().

    Resamples x and y with replacement independently; recomputes the mean
    difference on each resample; returns (observed_diff, ci_lower, ci_upper).
    """
    rng = rng if rng is not None else np.random.default_rng()
    nx, ny = len(x), len(y)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        xb = rng.choice(x, size=nx, replace=True)
        yb = rng.choice(y, size=ny, replace=True)
        diffs[b] = xb.mean() - yb.mean()
    lo = float(np.quantile(diffs, alpha / 2))
    hi = float(np.quantile(diffs, 1 - alpha / 2))
    return float(x.mean() - y.mean()), lo, hi


def cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    """Pooled Cohen's d. Returns 0 if pooled variance is 0."""
    nx, ny = len(x), len(y)
    if nx + ny <= 2:
        return 0.0
    vx, vy = x.var(ddof=1), y.var(ddof=1)
    pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    if pooled == 0:
        return 0.0
    return float((x.mean() - y.mean()) / pooled)


def permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = 1000,
    rng: np.random.Generator | None = None,
) -> float:
    """Two-sided permutation test on mean difference. Returns p-value."""
    rng = rng if rng is not None else np.random.default_rng()
    observed = abs(x.mean() - y.mean())
    combined = np.concatenate([np.asarray(x), np.asarray(y)])
    nx = len(x)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(combined)
        if abs(combined[:nx].mean() - combined[nx:].mean()) >= observed:
            count += 1
    return (count + 1) / (n_perm + 1)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR. Returns q-values aligned to input order.

    As bh_fdr in NP_MT_FALSIFY4_1_39.m, NaN (untested) entries are left NaN and
    excluded from the number of tests m.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p) & (p >= 0) & (p <= 1))
    if valid.size == 0:
        return out
    order = np.argsort(p[valid], kind="stable")
    ranked = p[valid][order]
    m = ranked.size
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[valid[order]] = np.clip(q, 0.0, 1.0)
    return out


def bootstrap_mean_ci(
    replicates: np.ndarray,
    n_boot: int = 10_000,
    alpha: float = 0.05,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Percentile bootstrap CI of the column means (bootstrap_mean_CI in MATLAB).

    ``replicates`` is (n_reps, n_nodes). Returns (ci_low, ci_high), each
    (n_nodes,). The intervals describe uncertainty in the estimated
    condition-specific mean across simulation initialisations, so they may be
    asymmetric for skewed (e.g. bistable) outcome distributions.
    """
    rng = rng if rng is not None else np.random.default_rng()
    x = np.asarray(replicates, dtype=float)
    n_nodes = x.shape[1]
    lo = np.full(n_nodes, np.nan)
    hi = np.full(n_nodes, np.nan)
    for j in range(n_nodes):
        col = x[:, j][np.isfinite(x[:, j])]
        if col.size == 0:
            continue
        if col.size == 1:
            lo[j] = hi[j] = col[0]
            continue
        draws = rng.integers(0, col.size, size=(n_boot, col.size))
        means = col[draws].mean(axis=1)
        lo[j], hi[j] = np.percentile(means, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return lo, hi


def _ranksum_safe(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided Wilcoxon rank-sum p (ranksum_safe in MATLAB)."""
    from scipy.stats import mannwhitneyu

    if a.size == 0 or b.size == 0:
        return float("nan")
    if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
        return 1.0
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def baseline_nonparametric_stats(
    hypo: np.ndarray,
    normal: np.ndarray,
    hyper: np.ndarray,
    node_names: list[str],
    exclude: tuple[str, ...] = ("Hypo", "NL", "HL"),
):
    """Kruskal–Wallis + pairwise rank-sum tests with BH-FDR, per node.

    Ports baseline_nonparametric_stats in NP_MT_FALSIFY4_1_39.m. Inputs are
    (n_reps, n_nodes) replicate steady states per regime. The fixed loading
    inputs are not tested. q-values are BH-adjusted across all tested nodes
    within each test family. These p/q values describe model-output variation
    among randomised initialisations, not biological replicates.

    Returns a pandas DataFrame with one row per node.
    """
    import pandas as pd
    from scipy.stats import kruskal

    groups = [np.asarray(g, dtype=float) for g in (hypo, normal, hyper)]
    n_nodes = len(node_names)
    p_kw = np.full(n_nodes, np.nan)
    p_hn = np.full(n_nodes, np.nan)
    p_hh = np.full(n_nodes, np.nan)
    p_nh = np.full(n_nodes, np.nan)
    med = np.full((n_nodes, 3), np.nan)
    q1 = np.full((n_nodes, 3), np.nan)
    q3 = np.full((n_nodes, 3), np.nan)
    count = np.zeros((n_nodes, 3), dtype=int)

    for j, name in enumerate(node_names):
        vals = [g[:, j][np.isfinite(g[:, j])] for g in groups]
        for c, v in enumerate(vals):
            count[j, c] = v.size
            if v.size:
                med[j, c] = np.median(v)
                q1[j, c], q3[j, c] = np.percentile(v, [25, 75])
        h, n, l = vals
        if name in exclude or min(h.size, n.size, l.size) == 0:
            continue
        pooled = np.concatenate([h, n, l])
        p_kw[j] = 1.0 if np.all(pooled == pooled[0]) else float(kruskal(h, n, l).pvalue)
        p_hn[j] = _ranksum_safe(h, n)
        p_hh[j] = _ranksum_safe(h, l)
        p_nh[j] = _ranksum_safe(n, l)

    return pd.DataFrame({
        "Node": list(node_names),
        "N_Hypo": count[:, 0], "N_NL": count[:, 1], "N_HL": count[:, 2],
        "Hypo_Median": med[:, 0], "Hypo_Q1": q1[:, 0], "Hypo_Q3": q3[:, 0],
        "NL_Median": med[:, 1], "NL_Q1": q1[:, 1], "NL_Q3": q3[:, 1],
        "HL_Median": med[:, 2], "HL_Q1": q1[:, 2], "HL_Q3": q3[:, 2],
        "p_KW_3group": p_kw, "q_KW_3group": bh_fdr(p_kw),
        "p_Hypo_NL": p_hn, "q_Hypo_NL": bh_fdr(p_hn),
        "p_Hypo_HL": p_hh, "q_Hypo_HL": bh_fdr(p_hh),
        "p_NL_HL": p_nh, "q_NL_HL": bh_fdr(p_nh),
    })


def permutation_pvalues(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = 1000,
    rng: np.random.RandomState | np.random.Generator | None = None,
) -> np.ndarray:
    """Column-wise two-sided permutation p-values for mean(x) - mean(y).

    Vectorised form of perm_p_two_sided in RESCUE_NEW4_1_48_targets_HQ.m, applied to
    every column (node) of x and y, shape (n_x, n_nodes) and (n_y, n_nodes).
    One set of n_perm label permutations is shared across columns; MATLAB draws
    a fresh randperm per node, so individual p-values differ at the level of
    permutation noise. p = (count + 1) / (n_perm + 1).
    """
    rng = rng if rng is not None else np.random.default_rng()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx = x.shape[0]
    pool = np.vstack([x, y])
    n = pool.shape[0]
    observed = np.abs(x.mean(axis=0) - y.mean(axis=0))
    count = np.zeros(pool.shape[1])
    for _ in range(n_perm):
        perm = rng.permutation(n)
        a = pool[perm[:nx]].mean(axis=0)
        b = pool[perm[nx:]].mean(axis=0)
        count += np.abs(a - b) >= observed
    return (count + 1) / (n_perm + 1)
