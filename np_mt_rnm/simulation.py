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

# The three chronic mechanical loading inputs. These are the only nodes that
# are externally driven rather than computed by the network, so they are the
# only SBML boundary species. Do NOT infer this set from "has no regulators":
# MATLAB gives any unregulated node omega = 0, so dX/dt = -X and it decays to
# 0 rather than being held (this bit NutD in the 353-edge 2026-09 revision).
MECHANICAL_INPUTS: tuple[str, ...] = ("Hypo", "NL", "HL")

# Per paper Section 2.2: "normal loading was represented by NL = 0.8
# (Hypo = HL = 0.01), hyper-loading by HL = 0.8 (Hypo = NL = 0.01), and
# hypo-loading by Hypo = 0.2 (NL = HL = 0.01)". Identical to
# RESCUE_NEW4_1_final.m:109-111. (NP_MT_RNM_FALSIFY4_1.m:53 uses NL = 0.10 for
# Hyper; we follow the paper. With seeded ensembles the falsification outcome
# is the same, 43/45, under either value.)
REGIME_PRESETS: dict[str, dict[str, float]] = {
    "Hypo":   {"Hypo": 0.20, "NL": 0.01, "HL": 0.01},
    "Normal": {"Hypo": 0.01, "NL": 0.80, "HL": 0.01},
    "Hyper":  {"Hypo": 0.01, "NL": 0.01, "HL": 0.80},
}

PAPER_REGIMES: tuple[str, ...] = ("Hypo", "Normal", "Hyper")

# Seed of replicate 1 for each baseline ensemble. RESCUE_NEW4_1_final.m seeds
# replicate r = 1..100 with rng(BASE_SEED + offset + r, 'twister') where
# BASE_SEED = 1 and offset = 1000 / 2000 / 3000 for Hypo / Normal / Hyper, so
# replicate r of e.g. Hyper uses seed 3001 + r. run_replicates gives replicate
# k (0-based) seed `seed + k`, so passing these values reproduces MATLAB's
# initial states exactly (see matlab_rand).
REGIME_SEEDS: dict[str, int] = {"Hypo": 1002, "Normal": 2002, "Hyper": 3002}

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


def matlab_rand(n: int, seed: int) -> np.ndarray:
    """Return MATLAB's ``rng(seed, 'twister'); rand(n, 1)`` for seed >= 1.

    MATLAB's default generator and numpy's legacy RandomState are the same
    Mersenne Twister (MT19937, init_genrand seeding, 53-bit doubles), so the
    draws are bit-identical. MATLAB maps seed 0 to 5489, so seeds below 1 do
    not correspond to MATLAB.
    """
    return np.random.RandomState(seed).random_sample(n)


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

    clamped_mask, x_clamp = build_clamps(net, regime=regime, user_clamps=user_clamps)

    if x0 is None:
        if seed is None:
            raise ValueError("either seed or x0 must be given")
        x0 = matlab_rand(n, seed)
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


def _one_replicate(net, regime, user_clamps, seed, x0=None):
    # Top-level helper so joblib can pickle it (when using loky backend).
    return run_single_replicate(
        net, regime=regime, user_clamps=user_clamps, seed=seed, x0=x0
    )


def run_replicates(
    net: Network,
    regime: Mapping[str, float],
    n_reps: int,
    user_clamps: Mapping[str, float] | None = None,
    seed: int = 0,
    n_jobs: int = -1,
    x0_list: np.ndarray | None = None,
) -> ReplicateEnsemble:
    """Run `n_reps` replicates in parallel and return the ensemble.

    Each replicate gets a deterministic seed = `seed + replicate_index`
    so the ensemble is fully reproducible regardless of `n_jobs`.

    `x0_list`, shape (n_reps, n_nodes), supplies an explicit initial state per
    replicate instead of drawing a random one. This is what makes the rescue
    screen paired-sequential: replicate r's perturbed run continues from
    replicate r's own Hyper steady state, mirroring
    run_reps_from_replicate_states_parallel in RESCUE_NEW4_1_final.m. When it
    is given, `seed` no longer influences the initial conditions.
    """
    seeds = [seed + k for k in range(n_reps)]

    if x0_list is None:
        starts: list[np.ndarray | None] = [None] * n_reps
    else:
        x0_arr = np.asarray(x0_list, dtype=float)
        if x0_arr.shape != (n_reps, len(net.node_names)):
            raise ValueError(
                f"x0_list must have shape ({n_reps}, {len(net.node_names)}), "
                f"got {x0_arr.shape}"
            )
        starts = [x0_arr[k] for k in range(n_reps)]

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_one_replicate)(net, regime, user_clamps, s, x0)
        for s, x0 in zip(seeds, starts)
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
