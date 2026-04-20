"""Hyper→Normal rescue screen. Ports legacy/RESCUE_NEW4_1.m.

Per paper Section 2.5:
  - Start from Hyper steady-state ensemble.
  - Apply single- or dual-node clamps (catabolic → 0, anabolic → 1).
  - Compute replicate-wise ΔX = X_perturbed − X_hyper_baseline.
  - Rank perturbations by mean |Δ| within a functional category.

Catabolic KD nodes: RhoA-E, PIEZO1, PI3K-E, FAK-E, ROS (clamp to 0).
Anabolic UP nodes: SOX9, PPARγ, HIF-1α, NRF2, IκBα (clamp to 1).

Total = 5 single KD + 5 single UP + 5×5 dual = 35 perturbations.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterator

import numpy as np

from np_mt_rnm.network import Network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates

# Using the Excel's exact Greek-letter node names.
CATABOLIC_DOWN_NODES: tuple[str, ...] = ("RhoA-E", "PIEZO1", "PI3K-E", "FAK-E", "ROS")
ANABOLIC_UP_NODES: tuple[str, ...] = ("SOX9", "PPARγ", "HIF-1α", "NRF2", "IκBα")


@dataclass(frozen=True)
class Perturbation:
    anabolic_up: str | None
    catabolic_down: str | None

    @property
    def label(self) -> str:
        parts: list[str] = []
        if self.anabolic_up is not None:
            parts.append(f"{self.anabolic_up}↑")
        if self.catabolic_down is not None:
            parts.append(f"{self.catabolic_down}↓")
        return " + ".join(parts) if parts else "baseline"


@dataclass(frozen=True)
class PerturbationResult:
    perturbation: Perturbation
    mean_delta: np.ndarray        # (n_nodes,) across replicates
    std_delta: np.ndarray         # (n_nodes,)
    node_names: list[str]
    n_reps: int


def enumerate_perturbations() -> Iterator[Perturbation]:
    for c in CATABOLIC_DOWN_NODES:
        yield Perturbation(anabolic_up=None, catabolic_down=c)
    for a in ANABOLIC_UP_NODES:
        yield Perturbation(anabolic_up=a, catabolic_down=None)
    for a, c in product(ANABOLIC_UP_NODES, CATABOLIC_DOWN_NODES):
        yield Perturbation(anabolic_up=a, catabolic_down=c)


def run_perturbation(
    net: Network,
    anabolic_up: str | None,
    catabolic_down: str | None,
    n_reps: int,
    seed: int = 0,
    n_jobs: int = -1,
) -> PerturbationResult:
    """Run a perturbation on top of the Hyper regime.

    Compares replicate-paired:
      - unperturbed Hyper baseline
      - Hyper + clamp(anabolic_up=1, catabolic_down=0)
    Returns the per-node mean delta and std delta across paired replicates.
    """
    hyper = REGIME_PRESETS["Hyper"]
    clamps: dict[str, float] = {}
    if anabolic_up is not None:
        clamps[anabolic_up] = 1.0
    if catabolic_down is not None:
        clamps[catabolic_down] = 0.0

    baseline = run_replicates(net, regime=hyper, n_reps=n_reps, seed=seed, n_jobs=n_jobs)
    perturbed = run_replicates(
        net, regime=hyper, n_reps=n_reps, seed=seed, n_jobs=n_jobs, user_clamps=clamps
    )
    delta = perturbed.steady_states - baseline.steady_states
    return PerturbationResult(
        perturbation=Perturbation(anabolic_up=anabolic_up, catabolic_down=catabolic_down),
        mean_delta=delta.mean(axis=0),
        std_delta=delta.std(axis=0, ddof=1),
        node_names=list(net.node_names),
        n_reps=n_reps,
    )


def mean_abs_displacement(result: PerturbationResult, nodes: list[str]) -> float:
    """Ranking metric: mean |Δ| over a set of nodes (a functional category)."""
    idx = [result.node_names.index(n) for n in nodes if n in result.node_names]
    if not idx:
        return 0.0
    return float(np.abs(result.mean_delta[idx]).mean())
