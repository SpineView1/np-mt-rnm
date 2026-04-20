"""Simulation runner — regimes, clamps, replicate loops.

Ports logic from legacy/NP_MT_RNM_FALSIFY4_1.m (sections that prepare
the loading regimes and execute ode45 per replicate).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.integrate import solve_ivp

from np_mt_rnm.network import Network
from np_mt_rnm.ode import squads_rhs

# Per paper Section 2.2 and MATLAB baseline values.
REGIME_PRESETS: dict[str, dict[str, float]] = {
    "Hypo":   {"Hypo": 0.20, "NL": 0.01, "HL": 0.01},
    "Normal": {"Hypo": 0.01, "NL": 0.80, "HL": 0.01},
    "Hyper":  {"Hypo": 0.01, "NL": 0.01, "HL": 0.80},
}

# Solver settings per spec Section 6 (match MATLAB odeset).
T_SPAN = (0.0, 100.0)
SOLVER_KWARGS = dict(
    method="RK45",
    rtol=1e-8,
    atol=1e-10,
    max_step=0.5,
)


def build_clamps(
    net: Network,
    regime: Mapping[str, float],
    user_clamps: Mapping[str, float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build (clamped_mask, x_clamp) arrays aligned to `net.node_names`.

    `regime` names (e.g. Hypo/NL/HL) and any `user_clamps` entries are
    marked clamped; their initial activation is set to the provided value
    and held there by the RHS.
    """
    n = len(net.node_names)
    clamped_mask = np.zeros(n, dtype=bool)
    x_clamp = np.zeros(n)
    name_to_idx = {name: i for i, name in enumerate(net.node_names)}

    for name, value in regime.items():
        if name not in name_to_idx:
            raise KeyError(f"regime input {name!r} not found in network nodes")
        idx = name_to_idx[name]
        clamped_mask[idx] = True
        x_clamp[idx] = float(value)

    if user_clamps:
        for name, value in user_clamps.items():
            if name not in name_to_idx:
                raise KeyError(f"clamp node {name!r} not found in network nodes")
            idx = name_to_idx[name]
            clamped_mask[idx] = True
            x_clamp[idx] = float(value)

    return clamped_mask, x_clamp


@dataclass(frozen=True)
class ReplicateResult:
    x_final: np.ndarray            # shape (n,)
    max_abs_derivative: float
    converged: bool
    t_final: float
    seed: int


def _random_initial_state(n: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform [0, 1] initial condition, matching MATLAB rand(n,1)."""
    return rng.uniform(0.0, 1.0, size=n)


def run_single_replicate(
    net: Network,
    regime: Mapping[str, float],
    user_clamps: Mapping[str, float] | None = None,
    seed: int | None = None,
    x0: np.ndarray | None = None,
) -> ReplicateResult:
    """Run a single ODE replicate to t=100 and return the steady state.

    Reproduces one iteration of the per-replicate loop in the MATLAB
    baseline script. Clamped nodes (regime inputs + user overrides) have
    their value enforced throughout the solve via the RHS mask and
    additionally overwritten in the final state for bit-exactness.
    """
    n = len(net.node_names)
    rng = np.random.default_rng(seed)

    clamped_mask, x_clamp = build_clamps(net, regime=regime, user_clamps=user_clamps)

    if x0 is None:
        x0 = _random_initial_state(n, rng)
    else:
        x0 = np.asarray(x0, dtype=float).copy()

    # Initialize clamped nodes exactly at their clamp values.
    x0[clamped_mask] = x_clamp[clamped_mask]

    def rhs(t, y):
        return squads_rhs(t, y, net.mact, net.minh, clamped=clamped_mask)

    sol = solve_ivp(rhs, T_SPAN, x0, **SOLVER_KWARGS)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")

    x_final = sol.y[:, -1].copy()
    # Clamped bit-exactness (guards against 1e-12 numerical drift).
    x_final[clamped_mask] = x_clamp[clamped_mask]

    dxdt_final = squads_rhs(sol.t[-1], x_final, net.mact, net.minh, clamped=clamped_mask)
    max_abs = float(np.abs(dxdt_final).max())

    return ReplicateResult(
        x_final=x_final,
        max_abs_derivative=max_abs,
        converged=max_abs < 1e-8,
        t_final=float(sol.t[-1]),
        seed=int(seed) if seed is not None else -1,
    )


from joblib import Parallel, delayed


@dataclass(frozen=True)
class ReplicateEnsemble:
    """Ensemble output of run_replicates."""

    steady_states: np.ndarray      # shape (n_reps, n_nodes)
    converged: np.ndarray          # shape (n_reps,), bool
    max_abs_derivatives: np.ndarray  # shape (n_reps,), float
    node_names: list[str]
    regime: dict[str, float]

    @property
    def all_converged(self) -> bool:
        return bool(self.converged.all())

    def mean(self) -> np.ndarray:
        return self.steady_states.mean(axis=0)

    def std(self) -> np.ndarray:
        return self.steady_states.std(axis=0, ddof=1)


def _one_replicate(net, regime, user_clamps, seed):
    # Top-level helper so joblib can pickle it (when using loky backend).
    return run_single_replicate(net, regime=regime, user_clamps=user_clamps, seed=seed)


def run_replicates(
    net: Network,
    regime: Mapping[str, float],
    n_reps: int,
    user_clamps: Mapping[str, float] | None = None,
    seed: int = 0,
    n_jobs: int = -1,
) -> ReplicateEnsemble:
    """Run `n_reps` replicates in parallel and return the ensemble.

    Each replicate gets a deterministic seed = `seed + replicate_index`
    so the ensemble is fully reproducible regardless of `n_jobs`.
    """
    seeds = [seed + k for k in range(n_reps)]

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_one_replicate)(net, regime, user_clamps, s) for s in seeds
    )

    steady = np.stack([r.x_final for r in results])
    conv = np.array([r.converged for r in results])
    maxd = np.array([r.max_abs_derivative for r in results])

    return ReplicateEnsemble(
        steady_states=steady,
        converged=conv,
        max_abs_derivatives=maxd,
        node_names=list(net.node_names),
        regime=dict(regime),
    )
