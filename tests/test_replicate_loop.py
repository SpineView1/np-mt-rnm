"""Tests for run_replicates — N-replicate ensemble at a given regime."""
from pathlib import Path

import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


@pytest.fixture(scope="module")
def net():
    return load_network(DATA_XLSX)


def test_run_replicates_returns_n_reps_rows(net):
    ens = run_replicates(net, regime=REGIME_PRESETS["Normal"], n_reps=4, seed=0, n_jobs=1)
    assert ens.steady_states.shape == (4, len(net.node_names))


def test_run_replicates_seed_reproducibility(net):
    ens1 = run_replicates(net, regime=REGIME_PRESETS["Normal"], n_reps=3, seed=42, n_jobs=1)
    ens2 = run_replicates(net, regime=REGIME_PRESETS["Normal"], n_reps=3, seed=42, n_jobs=1)
    np.testing.assert_allclose(ens1.steady_states, ens2.steady_states, atol=1e-10)


def test_run_replicates_all_converged(net):
    ens = run_replicates(net, regime=REGIME_PRESETS["Normal"], n_reps=4, seed=0, n_jobs=1)
    assert ens.all_converged, f"{(~ens.converged).sum()} replicates did not converge"
