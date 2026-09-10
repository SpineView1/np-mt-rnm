"""Tests for the TRUE Hyper->Normal rescue metric.

Ports analyze_group_true_rescue in legacy/RESCUE_NEW4_1_final.m:

    D_Hyper       = norm(hyper_profile - normal_profile)
    D_Rescue(c)   = norm(rescued_profile(:,c) - normal_profile)
    RescuePercent = 100 * (D_Hyper - D_Rescue) / D_Hyper

100% = the perturbation lands exactly on the Normal profile,
  0% = no better than untreated Hyper,
 <0% = the perturbation moves the group further from Normal.
"""
import numpy as np
import pytest

from np_mt_rnm.rescue import true_rescue_percent


def test_landing_on_normal_scores_100_percent():
    hyper = np.array([1.0, 1.0, 1.0])
    normal = np.array([0.0, 0.0, 0.0])
    rescued = normal.reshape(-1, 1)  # one strategy, exactly Normal
    pct = true_rescue_percent(hyper, normal, rescued)
    np.testing.assert_allclose(pct, [100.0])


def test_staying_at_hyper_scores_zero_percent():
    hyper = np.array([1.0, 1.0, 1.0])
    normal = np.array([0.0, 0.0, 0.0])
    rescued = hyper.reshape(-1, 1)
    pct = true_rescue_percent(hyper, normal, rescued)
    np.testing.assert_allclose(pct, [0.0], atol=1e-12)


def test_moving_away_from_normal_scores_negative():
    hyper = np.array([1.0, 0.0])
    normal = np.array([0.0, 0.0])
    rescued = np.array([[2.0], [0.0]])  # twice as far as Hyper
    pct = true_rescue_percent(hyper, normal, rescued)
    np.testing.assert_allclose(pct, [-100.0])


def test_halfway_scores_50_percent():
    hyper = np.array([2.0, 0.0])
    normal = np.array([0.0, 0.0])
    rescued = np.array([[1.0], [0.0]])
    pct = true_rescue_percent(hyper, normal, rescued)
    np.testing.assert_allclose(pct, [50.0])


def test_multiple_strategies_are_scored_columnwise():
    """rescued_profiles is (n_group_nodes, n_strategies), as FINAL(rowIdx,:) is."""
    hyper = np.array([1.0, 1.0])
    normal = np.array([0.0, 0.0])
    rescued = np.array(
        [[0.0, 1.0, 2.0],
         [0.0, 1.0, 2.0]]
    )
    pct = true_rescue_percent(hyper, normal, rescued)
    np.testing.assert_allclose(pct, [100.0, 0.0, -100.0], atol=1e-12)


def test_identical_hyper_and_normal_profiles_is_rejected():
    """MATLAB warns and skips when dHyper < eps; we refuse rather than divide by 0."""
    hyper = np.array([0.5, 0.5])
    normal = np.array([0.5, 0.5])
    rescued = np.array([[0.5], [0.5]])
    with pytest.raises(ValueError, match="indistinguishable"):
        true_rescue_percent(hyper, normal, rescued)


def test_metric_is_euclidean_not_per_node_mean():
    """Guards against accidentally using a mean-absolute distance instead.

    With a group whose deviation sits entirely in one node, the Euclidean and
    mean-absolute formulations give different answers; pin the Euclidean one.
    """
    hyper = np.array([3.0, 4.0])       # norm 5 from Normal
    normal = np.array([0.0, 0.0])
    rescued = np.array([[3.0], [0.0]])  # norm 3 from Normal
    pct = true_rescue_percent(hyper, normal, rescued)
    np.testing.assert_allclose(pct, [100.0 * (5.0 - 3.0) / 5.0])


def _tiny_inputs():
    node_names = ["A", "B", "C"]
    hyper = np.array([1.0, 1.0, 9.0])
    normal = np.array([0.0, 0.0, 9.0])
    # three strategies over group {A, B}: exact, no-op, worse
    final = np.array([
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        [9.0, 9.0, 9.0],
    ])
    return node_names, hyper, normal, final, ["exact", "noop", "worse"]


def test_rank_group_rescue_orders_best_first():
    from np_mt_rnm.rescue import rank_group_rescue

    names, hyper, normal, final, labels = _tiny_inputs()
    r = rank_group_rescue(
        "demo", ["A", "B"], names, hyper, normal, final, labels
    )
    assert r.strategy_labels == ["exact", "noop", "worse"]
    np.testing.assert_allclose(r.rescue_percent, [100.0, 0.0, -100.0], atol=1e-12)


def test_rank_group_rescue_respects_top_n():
    from np_mt_rnm.rescue import rank_group_rescue

    names, hyper, normal, final, labels = _tiny_inputs()
    r = rank_group_rescue(
        "demo", ["A", "B"], names, hyper, normal, final, labels, top_n=2
    )
    assert r.strategy_labels == ["exact", "noop"]
    assert r.rescue_percent.shape == (2,)


def test_rank_group_rescue_reports_missing_nodes():
    """Mirrors the MATLAB 'missing nodes' printout instead of silently dropping."""
    from np_mt_rnm.rescue import rank_group_rescue

    names, hyper, normal, final, labels = _tiny_inputs()
    r = rank_group_rescue(
        "demo", ["A", "B", "NoSuchNode"], names, hyper, normal, final, labels
    )
    assert r.used_nodes == ["A", "B"]
    assert r.missing_nodes == ["NoSuchNode"]


def test_rank_group_rescue_ignores_nodes_outside_the_group():
    """Node C sits at its Normal value in every strategy and must not dilute."""
    from np_mt_rnm.rescue import rank_group_rescue

    names, hyper, normal, final, labels = _tiny_inputs()
    r = rank_group_rescue("demo", ["A", "B"], names, hyper, normal, final, labels)
    assert r.used_nodes == ["A", "B"]
    np.testing.assert_allclose(r.distance_hyper_to_normal, np.sqrt(2.0))


def test_rank_group_rescue_rejects_label_count_mismatch():
    from np_mt_rnm.rescue import rank_group_rescue

    names, hyper, normal, final, _ = _tiny_inputs()
    with pytest.raises(ValueError, match="labels"):
        rank_group_rescue("demo", ["A", "B"], names, hyper, normal, final, ["only-one"])
