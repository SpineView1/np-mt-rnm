"""Constrained regime transition paths. Ports legacy/NP_MT_RNM_FSA4_1.m.

Two discrete paths per paper Section 2.3:
  Hypo → Normal: decrease Hypo 0.35→0.10, increase NL 0.10→0.35, HL fixed at 0.01
  Normal → Hyper: decrease NL  0.35→0.10, increase HL 0.10→0.35, Hypo fixed at 0.01

Each path has 6 discrete steps; at each step we run `n_reps` replicates and
record the mean activation per node.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from np_mt_rnm.network import Network
from np_mt_rnm.simulation import run_replicates

HYPO_TO_NORMAL_PATH: list[dict[str, float]] = [
    {"Hypo": 0.35, "NL": 0.10, "HL": 0.01},
    {"Hypo": 0.30, "NL": 0.15, "HL": 0.01},
    {"Hypo": 0.25, "NL": 0.20, "HL": 0.01},
    {"Hypo": 0.20, "NL": 0.25, "HL": 0.01},
    {"Hypo": 0.15, "NL": 0.30, "HL": 0.01},
    {"Hypo": 0.10, "NL": 0.35, "HL": 0.01},
]

NORMAL_TO_HYPER_PATH: list[dict[str, float]] = [
    {"Hypo": 0.01, "NL": 0.35, "HL": 0.10},
    {"Hypo": 0.01, "NL": 0.30, "HL": 0.15},
    {"Hypo": 0.01, "NL": 0.25, "HL": 0.20},
    {"Hypo": 0.01, "NL": 0.20, "HL": 0.25},
    {"Hypo": 0.01, "NL": 0.15, "HL": 0.30},
    {"Hypo": 0.01, "NL": 0.10, "HL": 0.35},
]


@dataclass(frozen=True)
class TransitionResult:
    path: list[dict[str, float]]
    means: np.ndarray                 # (n_steps, n_nodes)
    stds: np.ndarray
    node_names: list[str]


def run_transition_path(
    net: Network,
    path: Sequence[dict[str, float]],
    n_reps: int,
    seed: int = 0,
    n_jobs: int = -1,
) -> TransitionResult:
    """Run all steps in a path and return per-step mean/std activations."""
    means = np.zeros((len(path), len(net.node_names)))
    stds = np.zeros_like(means)
    for k, regime in enumerate(path):
        ens = run_replicates(
            net, regime=regime, n_reps=n_reps, seed=seed + 1000 * k, n_jobs=n_jobs
        )
        means[k] = ens.mean()
        stds[k] = ens.std()
    return TransitionResult(
        path=list(path),
        means=means,
        stds=stds,
        node_names=list(net.node_names),
    )
