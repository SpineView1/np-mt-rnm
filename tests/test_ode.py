"""Tests for np_mt_rnm.ode — SQUADS right-hand side.

Targets: np_mt_rnm/ode.py (ports legacy/ODESysFunS.m).

Expected values computed by hand in tests/ode_hand_computation.py.
Replace the two sentinels below with the actual printed values.
"""
import numpy as np
import pytest

from np_mt_rnm.ode import squads_rhs


# ⬇ Paste the actual values from `python tests/ode_hand_computation.py`.
_EXPECTED_A = -0.3298962834548918   # Case A (single activator)
_EXPECTED_B = 0.078217229825365109   # Case B (activator + inhibitor)


def _make_small_network(n_nodes: int, mact_edges: list[tuple[int, int]], minh_edges: list[tuple[int, int]]):
    mact = np.zeros((n_nodes, n_nodes), dtype=np.int8)
    minh = np.zeros((n_nodes, n_nodes), dtype=np.int8)
    for i, j in mact_edges:
        mact[i, j] = 1
    for i, j in minh_edges:
        minh[i, j] = 1
    return mact, minh


def test_squads_rhs_single_activator_case_A():
    """Case A: node 0 has single activator (node 1). X_1 = 0.6, X_0 = 0.4."""
    mact, minh = _make_small_network(2, mact_edges=[(0, 1)], minh_edges=[])
    x = np.array([0.4, 0.6])
    clamped = np.zeros(2, dtype=bool)
    dxdt = squads_rhs(0.0, x, mact, minh, clamped=clamped)
    assert dxdt.shape == (2,)
    np.testing.assert_allclose(dxdt[0], _EXPECTED_A, rtol=0, atol=1e-12)


def test_squads_rhs_activator_plus_inhibitor_case_B():
    """Case B: node 0 activated by node 1 and inhibited by node 2."""
    mact, minh = _make_small_network(3, mact_edges=[(0, 1)], minh_edges=[(0, 2)])
    x = np.array([0.3, 0.7, 0.2])
    clamped = np.zeros(3, dtype=bool)
    dxdt = squads_rhs(0.0, x, mact, minh, clamped=clamped)
    np.testing.assert_allclose(dxdt[0], _EXPECTED_B, rtol=0, atol=1e-12)


def test_squads_rhs_no_regulators_gives_only_decay():
    """A node with no activators or inhibitors has ω=0; dx/dt = RHS(ω=0) − γx.

    With the positive-exponent SQUADS formula used here, ω=0 gives activation=1,
    so dxdt = 1 - γ·X = 1 - 0.5 = 0.5. (Note: MATLAB ODESysFunS.m uses the
    negative-exponent variant where ω=0 → activation=0 → dxdt = -γX.)
    """
    mact = np.zeros((1, 1), dtype=np.int8)
    minh = np.zeros((1, 1), dtype=np.int8)
    x = np.array([0.5])
    clamped = np.zeros(1, dtype=bool)
    dxdt = squads_rhs(0.0, x, mact, minh, clamped=clamped)
    np.testing.assert_allclose(dxdt[0], 0.5, atol=1e-12)


def test_squads_rhs_clamped_nodes_have_zero_derivative():
    mact = np.zeros((2, 2), dtype=np.int8); mact[0, 1] = 1
    minh = np.zeros((2, 2), dtype=np.int8)
    x = np.array([0.4, 0.6])
    clamped = np.array([False, True])   # node 1 is clamped
    dxdt = squads_rhs(0.0, x, mact, minh, clamped=clamped)
    assert dxdt[1] == 0.0
