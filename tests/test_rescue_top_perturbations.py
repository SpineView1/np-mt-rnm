"""Encode paper Figs 10B/12B claims: top-ranked ECM rescuer = SOX9, top-ranked TF = NRF2."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_top_ecm_rescuer_involves_sox9():
    """Paper Fig 10B: the top-ranked ECM perturbation in the rescue screen
    involves SOX9 activation. The top-3 should contain at least one with SOX9."""
    df = pd.read_csv(ROOT / "results" / "tables" / "rescue_summary.csv")
    top = df.nlargest(3, "mean_abs_delta_ecm")
    assert any("SOX9" in p for p in top["perturbation"]), (
        f"expected SOX9 in top-3 ECM perturbations, got {list(top['perturbation'])}"
    )


def test_top_tf_rescuer_involves_nrf2():
    """Paper Fig 12B: the top-ranked TF perturbation involves NRF2 activation."""
    df = pd.read_csv(ROOT / "results" / "tables" / "rescue_summary.csv")
    top = df.nlargest(3, "mean_abs_delta_tf")
    assert any("NRF2" in p for p in top["perturbation"]), (
        f"expected NRF2 in top-3 TF perturbations, got {list(top['perturbation'])}"
    )
