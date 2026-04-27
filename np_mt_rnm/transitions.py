"""Constrained regime transition paths. Ports legacy/NP_MT_RNM_FSA4_1.m.

Two discrete paths per paper Section 2.3:
  Hypo → Normal: decrease Hypo 0.35→0.10, increase NL 0.10→0.35, HL fixed at 0.01
  Normal → Hyper: decrease NL  0.35→0.10, increase HL 0.10→0.35, Hypo fixed at 0.01

Each path has 6 discrete steps; at each step we run `n_reps` replicates and
record the mean activation per node.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from joblib import Parallel, delayed
from scipy.integrate import solve_ivp

from np_mt_rnm.network import Network
from np_mt_rnm.ode import squads_rhs
from np_mt_rnm.simulation import _random_initial_state, build_clamps, run_replicates

HYPO_TO_NORMAL_PATH: list[dict[str, float]] = [
    {"Hypo": 0.35, "NL": 0.10, "HL": 0.01},
    {"Hypo": 0.30, "NL": 0.15, "HL": 0.01},
    {"Hypo": 0.25, "NL": 0.20, "HL": 0.01},
    {"Hypo": 0.20, "NL": 0.25, "HL": 0.01},
    {"Hypo": 0.15, "NL": 0.30, "HL": 0.01},
    {"Hypo": 0.10, "NL": 0.35, "HL": 0.01},
]

NORMAL_TO_HYPER_PATH: list[dict[str, float]] = [
    {"Hypo": 0.01, "NL": 0.35, "HL": 0.10},
    {"Hypo": 0.01, "NL": 0.30, "HL": 0.15},
    {"Hypo": 0.01, "NL": 0.25, "HL": 0.20},
    {"Hypo": 0.01, "NL": 0.20, "HL": 0.25},
    {"Hypo": 0.01, "NL": 0.15, "HL": 0.30},
    {"Hypo": 0.01, "NL": 0.10, "HL": 0.35},
]


@dataclass(frozen=True)
class TransitionResult:
    path: list[dict[str, float]]
    means: np.ndarray                 # (n_steps, n_nodes)
    stds: np.ndarray
    node_names: list[str]


def run_transition_path(
    net: Network,
    path: Sequence[dict[str, float]],
    n_reps: int,
    seed: int = 0,
    n_jobs: int = -1,
) -> TransitionResult:
    """Run all steps in a path and return per-step mean/std activations."""
    means = np.zeros((len(path), len(net.node_names)))
    stds = np.zeros_like(means)
    for k, regime in enumerate(path):
        ens = run_replicates(
            net, regime=regime, n_reps=n_reps, seed=seed + 1000 * k, n_jobs=n_jobs
        )
        means[k] = ens.mean()
        stds[k] = ens.std()
    return TransitionResult(
        path=list(path),
        means=means,
        stds=stds,
        node_names=list(net.node_names),
    )


@dataclass(frozen=True)
class TrajectoryResult:
    t: np.ndarray                 # (n_t,)
    mean: np.ndarray              # (n_t, n_nodes)
    std: np.ndarray               # (n_t, n_nodes)
    node_names: list[str]
    n_reps_used: int


def _one_trajectory(
    net: Network,
    regime_a: Mapping[str, float],
    regime_b: Mapping[str, float],
    t_switch: float,
    t_end: float,
    t_eval: np.ndarray,
    seed: int,
) -> np.ndarray:
    """Single replicate trajectory across the regime switch. Returns (n_t, n_nodes)."""
    rng = np.random.default_rng(seed)
    n = len(net.node_names)

    mask_a, x_clamp_a = build_clamps(net, regime=regime_a)
    mask_b, x_clamp_b = build_clamps(net, regime=regime_b)

    x0 = _random_initial_state(n, rng)
    x0[mask_a] = x_clamp_a[mask_a]

    t_eval_a = t_eval[t_eval <= t_switch]
    t_eval_b = t_eval[t_eval > t_switch]

    def rhs_a(t, y):
        return squads_rhs(t, y, net.mact, net.minh, clamped=mask_a)

    def rhs_b(t, y):
        return squads_rhs(t, y, net.mact, net.minh, clamped=mask_b)

    sol_a = solve_ivp(
        rhs_a, (0.0, t_switch), x0,
        t_eval=t_eval_a, method="LSODA",
        rtol=1e-6, atol=1e-9,
    )
    if not sol_a.success:
        raise RuntimeError(f"phase A solve failed: {sol_a.message}")
    y_at_switch = sol_a.y[:, -1].copy()
    y_at_switch[mask_b] = x_clamp_b[mask_b]

    sol_b = solve_ivp(
        rhs_b, (t_switch, t_end), y_at_switch,
        t_eval=t_eval_b, method="LSODA",
        rtol=1e-6, atol=1e-9,
    )
    if not sol_b.success:
        raise RuntimeError(f"phase B solve failed: {sol_b.message}")

    out = np.zeros((len(t_eval), n))
    out[: len(t_eval_a)] = sol_a.y.T
    out[len(t_eval_a) :] = sol_b.y.T
    return out


def run_trajectory(
    net: Network,
    regime_a: Mapping[str, float],
    regime_b: Mapping[str, float],
    t_switch: float,
    t_end: float,
    n_t: int,
    n_reps: int,
    seed: int = 0,
    n_jobs: int = 1,
) -> TrajectoryResult:
    """Run replicate trajectories with a regime switch at t_switch.

    Phase A: solve under regime_a from t=0 to t_switch.
    Phase B: continue under regime_b from the switch state to t_end.
    Both phases use t_eval from a uniform grid of n_t samples in [0, t_end].
    """
    if t_switch <= 0 or t_switch >= t_end:
        raise ValueError("require 0 < t_switch < t_end")
    t_eval = np.linspace(0.0, t_end, n_t)
    runs = Parallel(n_jobs=n_jobs)(
        delayed(_one_trajectory)(
            net, regime_a, regime_b, t_switch, t_end, t_eval, seed + 7919 * k
        )
        for k in range(n_reps)
    )
    stacked = np.stack(runs, axis=0)
    return TrajectoryResult(
        t=t_eval,
        mean=stacked.mean(axis=0),
        std=stacked.std(axis=0, ddof=1) if n_reps > 1 else np.zeros_like(stacked[0]),
        node_names=list(net.node_names),
        n_reps_used=n_reps,
    )
