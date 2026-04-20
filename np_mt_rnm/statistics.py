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
    """Benjamini–Hochberg FDR. Returns q-values aligned to input order."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    out = np.empty(n)
    out[order] = np.clip(q, 0.0, 1.0)
    return out
