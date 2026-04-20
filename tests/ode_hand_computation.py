"""Hand-computed SQUADS reference values used in tests/test_ode.py.

Run this file directly to see the expected dX/dt values:
    python tests/ode_hand_computation.py
"""
import math

H = 10.0
GAMMA = 1.0
EXP_HALF_H = math.exp(5.0)  # exp(0.5 * h)


def _activation_term(omega: float) -> float:
    """The SQUADS activation function value at a given omega.

    Matches MATLAB ODESysFunS.m line 77 (negative exponent).
    """
    exp_h_shift = math.exp(-H * (omega - 0.5))
    return (-EXP_HALF_H + exp_h_shift) / ((1.0 - EXP_HALF_H) * (1.0 + exp_h_shift))


# --- Case A: node 0 has one activator (node 1), no inhibitors.
# α = 1 for the single activator, X_1 = 0.6, X_0 = 0.4.
#   A = ((1 + α)/α) * (α * X_1 / (1 + α * X_1))
#     = (2/1) * (0.6 / 1.6) = 0.75
# ω = A (only activators present)
A_case = (1 + 1) / 1 * (1 * 0.6) / (1 + 1 * 0.6)
omega_A = A_case
dxdt_A = _activation_term(omega_A) - GAMMA * 0.4
print(f"Case A (single activator): ω = {omega_A:.17g}, dX_0/dt = {dxdt_A:.17g}")

# --- Case B: node 0 has one activator (node 1, X_1 = 0.7) AND one inhibitor (node 2, X_2 = 0.2).
# X_0 = 0.3 (current activation of target node).
#   A = (2/1) * (0.7 / 1.7)
#   I_inner = (2/1) * (0.2 / 1.2)
#   I_term  = 1 - I_inner
#   ω = A * I_term
A_B = (1 + 1) / 1 * (1 * 0.7) / (1 + 1 * 0.7)
I_inner_B = (1 + 1) / 1 * (1 * 0.2) / (1 + 1 * 0.2)
I_term_B = 1.0 - I_inner_B
omega_B = A_B * I_term_B
dxdt_B = _activation_term(omega_B) - GAMMA * 0.3
print(f"Case B (activator+inhibitor): ω = {omega_B:.17g}, dX_0/dt = {dxdt_B:.17g}")
