"""Tests for simulation.py helpers: regime presets, clamp construction."""
import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import (
    REGIME_PRESETS,
    build_clamps,
)
from pathlib import Path

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


def test_regime_presets_match_matlab():
    """Mirrors the clamp constants in the two MATLAB scripts.

    NP_MT_RNM_FALSIFY4_1.m:51-53 sets the baseline/falsification regimes and
    RESCUE_NEW4_1_final.m:109-111 sets the rescue ones. They agree on Hypo and
    Normal but differ on NL in Hyper (0.10 vs 0.01), so both are kept.
    """
    assert REGIME_PRESETS["Hypo"] == {"Hypo": 0.20, "NL": 0.01, "HL": 0.01}
    assert REGIME_PRESETS["Normal"] == {"Hypo": 0.01, "NL": 0.80, "HL": 0.01}
    assert REGIME_PRESETS["Hyper"] == {"Hypo": 0.01, "NL": 0.10, "HL": 0.80}
    assert REGIME_PRESETS["Hyper_rescue"] == {"Hypo": 0.01, "NL": 0.01, "HL": 0.80}


def test_build_clamps_marks_regime_inputs():
    net = load_network(DATA_XLSX)
    clamped_mask, x_clamp = build_clamps(
        net, regime=REGIME_PRESETS["Normal"], user_clamps=None
    )
    for name in ("Hypo", "NL", "HL"):
        idx = net.node_names.index(name)
        assert clamped_mask[idx], f"{name} should be clamped"
    nl_idx = net.node_names.index("NL")
    assert x_clamp[nl_idx] == 0.80


def test_build_clamps_honors_user_clamps():
    net = load_network(DATA_XLSX)
    clamped_mask, x_clamp = build_clamps(
        net, regime=REGIME_PRESETS["Hyper"], user_clamps={"SOX9": 1.0}
    )
    sox9_idx = net.node_names.index("SOX9")
    assert clamped_mask[sox9_idx]
    assert x_clamp[sox9_idx] == 1.0


def test_build_clamps_unknown_user_node_raises():
    net = load_network(DATA_XLSX)
    with pytest.raises(KeyError):
        build_clamps(net, regime=REGIME_PRESETS["Normal"], user_clamps={"NOT_A_NODE": 1})


def test_run_replicates_honors_supplied_initial_states():
    """`x0_list` must seed each replicate, not be ignored in favour of RNG.

    Required by the MATLAB paired-sequential rescue design, where replicate r's
    perturbed run starts from replicate r's own Hyper steady state
    (run_reps_from_replicate_states_parallel in RESCUE_NEW4_1_final.m).
    """
    import numpy as np
    from np_mt_rnm.network import load_network
    from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

    net = load_network(DATA_XLSX)
    n = len(net.node_names)

    first = run_replicates(net, regime=REGIME_PRESETS["Hyper"], n_reps=3, seed=0, n_jobs=1)
    # Re-integrating a converged state must leave it where it is.
    second = run_replicates(
        net,
        regime=REGIME_PRESETS["Hyper"],
        n_reps=3,
        seed=999,                      # different seed: ignored when x0_list given
        n_jobs=1,
        x0_list=first.steady_states,
    )
    assert second.steady_states.shape == (3, n)
    np.testing.assert_allclose(second.steady_states, first.steady_states, atol=1e-6)
