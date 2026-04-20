"""SQUADS ODE right-hand side. Ports legacy/ODESysFunS.m lines 37–99.

The SQUADS formalism (Standardized Qualitative Dynamical Systems) is:

    dX_i/dt = (-exp(0.5h) + exp(h(ω_i - 0.5)))
              / ((1 - exp(0.5h)) * (1 + exp(h(ω_i - 0.5))))
              - γ_i * X_i

where ω_i ∈ [0, 1] is the effective regulatory input for node i,
combining activator and inhibitor contributions:

    A_i = ((1 + ΣA)/ΣA) * (ΣA·X_active)/(1 + ΣA·X_active)      (activators)
    I_i = ((1 + ΣB)/ΣB) * (ΣB·X_inhib)/(1 + ΣB·X_inhib)        (inhibitors)

    ω_i = A_i · (1 - I_i)    if both activators and inhibitors present
          A_i                if only activators present
          1 - I_i            if only inhibitors present
          0                  if neither present

All weights α, β are uniform = 1 in this model. γ = 1, h = 10.

Clamped nodes (Hypo, NL, HL, any user-specified overrides) have their
derivative forced to 0, matching MATLAB wrapper ODESysFun_Clamped_safe3.
"""
from __future__ import annotations

import numpy as np

H = 10.0
GAMMA = 1.0
# Precompute constants used for every RHS call.
_EXP_HALF_H = np.exp(0.5 * H)   # exp(5)
_DENOM_ONE_MINUS = 1.0 - _EXP_HALF_H    # (1 - exp(5))


def squads_rhs(
    t: float,
    x: np.ndarray,
    mact: np.ndarray,
    minh: np.ndarray,
    clamped: np.ndarray | None = None,
) -> np.ndarray:
    """Return dX/dt for the SQUADS system.

    Parameters
    ----------
    t : float
        Time (unused — system is autonomous). Accepted for solve_ivp signature.
    x : (n,) float array
        Current activation vector, values in [0, 1].
    mact, minh : (n, n) int arrays
        Signed binary adjacency matrices. mact[i, j] = 1 means j activates i.
    clamped : (n,) bool array, optional
        If True for index i, dX_i/dt is forced to 0.

    Returns
    -------
    (n,) float array of derivatives.
    """
    n = x.shape[0]
    if clamped is None:
        clamped = np.zeros(n, dtype=bool)

    # Activator and inhibitor contributions per target node.
    # A[i] = ((1 + k_A)/k_A) * (sum_j α·X_j) / (1 + sum_j α·X_j)
    #      where k_A = number of activators of i; α = 1.
    sum_act_x = mact @ x        # Σ α X over activators of each target (α=1)
    k_act = mact.sum(axis=1)    # number of activators per target
    sum_inh_x = minh @ x
    k_inh = minh.sum(axis=1)

    # Compute A_i where it applies (k_act > 0); zero otherwise.
    with np.errstate(divide="ignore", invalid="ignore"):
        a_scale = np.where(k_act > 0, (1.0 + k_act) / np.where(k_act > 0, k_act, 1), 0.0)
        a_term = np.where(
            k_act > 0,
            a_scale * (sum_act_x / (1.0 + sum_act_x)),
            0.0,
        )
        i_scale = np.where(k_inh > 0, (1.0 + k_inh) / np.where(k_inh > 0, k_inh, 1), 0.0)
        i_inner = np.where(
            k_inh > 0,
            i_scale * (sum_inh_x / (1.0 + sum_inh_x)),
            0.0,
        )

    # ω per regulation regime (vectorized via boolean masks).
    has_act = k_act > 0
    has_inh = k_inh > 0
    both = has_act & has_inh
    only_act = has_act & ~has_inh
    only_inh = ~has_act & has_inh

    omega = np.zeros(n)
    omega[both] = a_term[both] * (1.0 - i_inner[both])
    omega[only_act] = a_term[only_act]
    omega[only_inh] = 1.0 - i_inner[only_inh]
    # nodes with neither: omega stays 0

    # SQUADS activation curve.
    exp_h_shift = np.exp(H * (omega - 0.5))
    numer = -_EXP_HALF_H + exp_h_shift
    denom = _DENOM_ONE_MINUS * (1.0 + exp_h_shift)
    activation = numer / denom

    dxdt = activation - GAMMA * x
    dxdt[clamped] = 0.0
    return dxdt
