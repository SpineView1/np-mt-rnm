"""Baseline statistics helpers ported from NP_MT_FALSIFY4_1_39.m."""
import numpy as np

from np_mt_rnm.statistics import baseline_nonparametric_stats, bh_fdr, bootstrap_mean_ci


def test_bh_fdr_matches_hand_computation_and_skips_nan():
    p = np.array([0.01, np.nan, 0.04, 0.03, 0.20])
    q = bh_fdr(p)
    assert np.isnan(q[1])
    # m = 4 valid tests: sorted 0.01, 0.03, 0.04, 0.20 -> 0.04, 0.0533, 0.0533, 0.20
    np.testing.assert_allclose(q[[0, 3, 2, 4]], [0.04, 0.16 / 3, 0.16 / 3, 0.20])


def test_bootstrap_mean_ci_brackets_mean_and_collapses_for_constants():
    rng = np.random.default_rng(0)
    x = np.column_stack([rng.normal(0.5, 0.1, 100), np.full(100, 0.3)])
    lo, hi = bootstrap_mean_ci(x, n_boot=2000, rng=np.random.default_rng(1))
    assert lo[0] < x[:, 0].mean() < hi[0]
    assert hi[0] - lo[0] < 0.06
    assert lo[1] == hi[1] == 0.3


def test_nonparametric_stats_excludes_inputs_and_detects_shift():
    rng = np.random.default_rng(0)
    names = ["Hypo", "A", "B"]
    h = np.column_stack([np.full(30, 0.2), rng.random(30), np.full(30, 0.5)])
    n = np.column_stack([np.full(30, 0.01), rng.random(30) + 1, np.full(30, 0.5)])
    l = np.column_stack([np.full(30, 0.01), rng.random(30) + 2, np.full(30, 0.5)])
    st = baseline_nonparametric_stats(h, n, l, names, exclude=("Hypo",))
    assert np.isnan(st.loc[0, "p_KW_3group"])
    assert st.loc[1, "q_KW_3group"] < 1e-6 and st.loc[1, "q_NL_HL"] < 1e-6
    assert st.loc[2, "p_KW_3group"] == 1.0 and st.loc[2, "p_NL_HL"] == 1.0
