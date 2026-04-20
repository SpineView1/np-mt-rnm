"""Encode paper Fig 4–5 regime-polarity claims as test invariants.

These match what the biologist verifies by eye when running the MATLAB
scripts. Encoding them as tests lets CI catch regressions after the port.
"""
from pathlib import Path

import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


@pytest.fixture(scope="module")
def ensembles_and_net():
    net = load_network(DATA_XLSX)
    ens = {
        regime: run_replicates(net, regime=preset, n_reps=25, seed=1000, n_jobs=-1)
        for regime, preset in REGIME_PRESETS.items()
    }
    return ens, net


def _mean(net, ensemble, node):
    idx = net.node_names.index(node)
    return float(ensemble.mean()[idx])


def test_acan_high_under_normal_low_under_hyper(ensembles_and_net):
    ens, net = ensembles_and_net
    assert _mean(net, ens["Normal"], "ACAN") > 0.8
    assert _mean(net, ens["Hyper"], "ACAN") < 0.3


def test_col10a1_low_under_normal_high_under_hyper(ensembles_and_net):
    ens, net = ensembles_and_net
    assert _mean(net, ens["Normal"], "COL10A1") < 0.1
    assert _mean(net, ens["Hyper"], "COL10A1") > 0.8


def test_mmp13_low_under_normal_high_under_hyper(ensembles_and_net):
    ens, net = ensembles_and_net
    assert _mean(net, ens["Normal"], "MMP13") < 0.1
    assert _mean(net, ens["Hyper"], "MMP13") > 0.8


def test_inflammatory_low_under_normal(ensembles_and_net):
    """TNF and IL6 should stay near zero under Normal loading."""
    ens, net = ensembles_and_net
    for cyt in ("TNF", "IL6"):
        assert cyt in net.node_names, f"{cyt} not in network — check Excel"
        assert _mean(net, ens["Normal"], cyt) < 0.15, f"{cyt} too active under Normal"


def test_sox9_high_under_normal(ensembles_and_net):
    ens, net = ensembles_and_net
    assert _mean(net, ens["Normal"], "SOX9") > 0.7
