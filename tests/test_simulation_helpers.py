"""Tests for simulation.py helpers: regime presets, clamp construction."""
import numpy as np
import pytest

from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import (
    REGIME_PRESETS,
    REGIME_SEEDS,
    build_clamps,
    matlab_rand,
)
from pathlib import Path

DATA_XLSX = Path(__file__).resolve().parents[1] / "data" / "MT_PRIMARY4_1.xlsx"


def test_regime_presets_match_paper():
    """Paper Section 2.2 and RESCUE_NEW4_1_final.m:109-111.

    NP_MT_RNM_FALSIFY4_1.m:53 uses NL = 0.10 for Hyper; the paper states
    Hypo = NL = 0.01, which is what we follow.
    """
    assert REGIME_PRESETS == {
        "Hypo": {"Hypo": 0.20, "NL": 0.01, "HL": 0.01},
        "Normal": {"Hypo": 0.01, "NL": 0.80, "HL": 0.01},
        "Hyper": {"Hypo": 0.01, "NL": 0.01, "HL": 0.80},
    }


def test_regime_seeds_match_rescue_script():
    """rng(BASE_SEED + offset + r, 'twister') with BASE_SEED = 1, r = 1..100."""
    assert REGIME_SEEDS == {"Hypo": 1002, "Normal": 2002, "Hyper": 3002}


def test_matlab_rand_matches_matlab_twister():
    """First MT19937 doubles for init_genrand(1), i.e. rng(1,'twister'); rand(3,1).

    The end-to-end proof that the streams agree is tests/test_paper_reproduction.py:
    the rescue percentages only match the manuscript with these draws.
    """
    import numpy as np

    np.testing.assert_array_equal(
        matlab_rand(3, 1),
        [0.417022004702574, 0.7203244934421581, 0.00011437481734488664],
    )


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

    first = run_replicates(net, regime=REGIME_PRESETS["Hyper"], n_reps=3, seed=3002, n_jobs=1)
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
