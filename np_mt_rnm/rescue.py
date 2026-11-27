"""Hyper→Normal rescue screen. Ports legacy/RESCUE_NEW4_1_48_targets_HQ.m.

Per the manuscript (Methods, "Hyper-to-Normal rescue screening via node clamping"):
  - Start from the Hyper (HL) steady-state ensemble.
  - Apply single- or dual-node clamps (suppression → 0, activation → 1).
  - Compute replicate-wise ΔX = X_perturbed − X_hyper_baseline.
  - Score each functional module by TRUE rescue: the percentage reduction of
    the Euclidean HL-to-NL distance of the module's ensemble-mean profile.

Suppression targets: RhoA-E, PIEZO1, ROS, P38, mTORC1, NF-κB (clamp to 0).
Activation targets:  SOX9, PPARγ, HIF-1α, NRF2, SIRT1, FOXO (clamp to 1).

Total = 6 single suppressions + 6 single activations + 6×6 duals = 48.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterator

import numpy as np

from np_mt_rnm.network import Network
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, ReplicateEnsemble, run_replicates

# CATAB_KD_LIST / ANAB_UP_LIST in RESCUE_NEW4_1_48_targets_HQ.m, in its order,
# using the Excel's exact node names (P38 is uppercase in MT_PRIMARY4_1.xlsx).
CATABOLIC_DOWN_NODES: tuple[str, ...] = ("RhoA-E", "PIEZO1", "ROS", "P38", "mTORC1", "NF-κB")
ANABOLIC_UP_NODES: tuple[str, ...] = ("SOX9", "PPARγ", "HIF-1α", "NRF2", "SIRT1", "FOXO")

N_SINGLES = len(CATABOLIC_DOWN_NODES) + len(ANABOLIC_UP_NODES)          # 12
N_DUALS = len(CATABOLIC_DOWN_NODES) * len(ANABOLIC_UP_NODES)            # 36
N_PERTURBATIONS = N_SINGLES + N_DUALS                                   # 48

# Readout groups of the rescue script (its own group_categories, which differ
# slightly from the falsification script's display groups: e.g. the metabolic
# group omits the antioxidant enzymes and the ion-channel group omits AQP1/5).
# analyze_group_true_rescue resolves them with case-sensitive ismember, so a
# name that is not a network node (ATF2, NutD) is reported missing and skipped.
# Keys are the category ids of np_mt_rnm.categories.
RESCUE_GROUPS: dict[str, tuple[str, ...]] = {
    "cytokines_chemokines_proteases": (
        "TNF", "IL6", "IL1β", "IL8", "CCL2", "CXCL1", "CXCL3",
        "ADAMTS4/5", "MMP1", "MMP3", "MMP13",
    ),
    "growth_factor": (
        "TGFβ", "VEGF", "IGF1", "BMP2", "CCN2", "GDF5", "FGF2", "FGF18", "Wnt3a", "Wnt5a",
    ),
    "transcription_factor": (
        "CREB", "HIF-1α", "HIF-2α", "NF-κB", "AP-1", "FOXO", "SOX9", "ATF2", "NFAT", "RUNX2",
        "YAP/TAZ", "MRTF-A", "NRF2", "HSF1", "TonEBP", "ELK1", "PPARγ", "CITED2",
    ),
    "ion_channel": (
        "Ca2+os", "Ca2+su", "CaMKII", "PKC-E", "PKC-M", "PLCγ-M", "PLCγ-E", "CaN", "IP3", "PLA2",
    ),
    "metabolic": (
        "LKB1", "NAD+", "AMPK", "mTORC1", "mTORC2", "SIRT1", "PI3K-M", "PI3K-E", "PIP3-M",
        "PIP3-E", "PDK1-M", "PDK1-E", "AKT1-M", "AKT1-E", "GSK3B", "ULK1", "PTEN", "PLD2",
        "PGE2", "COX-2", "PHD2", "VHL", "Rheb",
    ),
    "ecm_matrix": ("COL2A1", "COL1A1", "COL10A1", "ACAN", "TIMP3"),
    "mechanosensor": (
        "Hypo", "NL", "HL", "α5β1-FN", "α5β1-Fs", "αVβ3", "αVβ6", "SDC4-E", "SDC4-M",
        "TRPV4", "PIEZO1", "NutD", "MitD",
    ),
    "cell_fate": ("Bcl2", "BAX", "CASP3", "CASP9", "BNIP3", "GADD45", "DRP1", "MOMP"),
    "oxidative_proteostasis": ("HO-1", "GPX1", "SOD1", "SOD2", "CAT", "HSP70", "HSP27", "ROS"),
    "mapk": (
        "RAS-M", "RAS-E", "RAF-M", "RAF-E", "MEK-M", "MEK-E", "ERK-M", "ERK-E",
        "MKK3/6", "MKK4/7", "JNK", "P38", "RSK", "TAK1",
    ),
    "rho_cytoskeletal": (
        "RhoA-M", "RhoA-E", "RAC1-M", "RAC1-E", "CDC42", "ROCK-M", "ROCK-E", "PAK1", "PKN1",
        "FAK-M", "FAK-E", "MST1/2", "LATS1/2",
    ),
}

# The eight modules whose Top-20 rankings the manuscript reports
# (main-text rescue figures and Supplementary Table "Highest-scoring
# intervention for different functional modules").
PAPER_RESCUE_MODULES: tuple[str, ...] = (
    "ecm_matrix",
    "growth_factor",
    "transcription_factor",
    "cytokines_chemokines_proteases",
    "metabolic",
    "ion_channel",
    "oxidative_proteostasis",
    "cell_fate",
)

# Phenotype readout panels (ANABOLIC_TX / CATABOLIC_TX in the rescue script).
ANABOLIC_PANEL: tuple[str, ...] = (
    "ACAN", "COL2A1", "TIMP3", "PPARγ", "GDF5", "IGF1", "SOX9", "Bcl2",
)
CATABOLIC_PANEL: tuple[str, ...] = (
    "MMP1", "MMP3", "MMP13", "ADAMTS4/5", "IL1β", "IL6", "IL8", "TNF", "iNOS", "COX-2",
)


@dataclass(frozen=True)
class Perturbation:
    anabolic_up: str | None
    catabolic_down: str | None

    @property
    def short_label(self) -> str:
        """MATLAB ShortClampLabels style, as printed in the manuscript: 'SOX9-ROS'."""
        parts = [n for n in (self.anabolic_up, self.catabolic_down) if n is not None]
        return "-".join(parts) if parts else "baseline"

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
    baseline_states: np.ndarray   # (n_reps, n_nodes) Hyper steady states
    perturbed_states: np.ndarray  # (n_reps, n_nodes) continued from baseline

    @property
    def baseline_mean(self) -> np.ndarray:
        return self.baseline_states.mean(axis=0)

    @property
    def perturbed_mean(self) -> np.ndarray:
        """Column of FINAL in RESCUE_NEW4_1_48_targets_HQ.m."""
        return self.perturbed_states.mean(axis=0)


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
    seed: int = REGIME_SEEDS["Hyper"],
    n_jobs: int = -1,
    baseline: ReplicateEnsemble | None = None,
) -> PerturbationResult:
    """Run a perturbation on top of the Hyper regime, paired-sequentially.

    Mirrors the design stated in RESCUE_NEW4_1_48_targets_HQ.m; for each replicate r:

        random initial state -> Hyper steady state -> perturbed steady state

    So the perturbed run *continues from* replicate r's own Hyper steady state
    rather than restarting from a fresh random state with the clamp applied at
    t=0. The delta is then replicate-wise, Delta_r = Perturbed_r - Hyper_r.

    `baseline` supplies an already-computed Hyper ensemble. MATLAB solves
    SS_hyper_all once and reuses it for all 48 perturbations, so the whole
    screen shares one baseline; pass it in to reproduce that (and halve the
    solver work). When omitted, a fresh baseline is solved for this call.
    """
    hyper = REGIME_PRESETS["Hyper"]
    clamps: dict[str, float] = {}
    if anabolic_up is not None:
        clamps[anabolic_up] = 1.0
    if catabolic_down is not None:
        clamps[catabolic_down] = 0.0

    if baseline is None:
        baseline = run_replicates(
            net, regime=hyper, n_reps=n_reps, seed=seed, n_jobs=n_jobs
        )
    elif baseline.steady_states.shape[0] != n_reps:
        raise ValueError(
            f"baseline has {baseline.steady_states.shape[0]} replicates, "
            f"n_reps={n_reps}"
        )
    perturbed = run_replicates(
        net,
        regime=hyper,
        n_reps=n_reps,
        seed=seed,
        n_jobs=n_jobs,
        user_clamps=clamps,
        x0_list=baseline.steady_states,
    )
    delta = perturbed.steady_states - baseline.steady_states
    return PerturbationResult(
        perturbation=Perturbation(anabolic_up=anabolic_up, catabolic_down=catabolic_down),
        mean_delta=delta.mean(axis=0),
        std_delta=delta.std(axis=0, ddof=1),
        node_names=list(net.node_names),
        n_reps=n_reps,
        baseline_states=baseline.steady_states,
        perturbed_states=perturbed.steady_states,
    )


def mean_abs_displacement(result: PerturbationResult, nodes: list[str]) -> float:
    """Ranking metric: mean |Δ| over a set of nodes (a functional category)."""
    idx = [result.node_names.index(n) for n in nodes if n in result.node_names]
    if not idx:
        return 0.0
    return float(np.abs(result.mean_delta[idx]).mean())


def true_rescue_percent(
    hyper_profile: np.ndarray,
    normal_profile: np.ndarray,
    rescued_profiles: np.ndarray,
) -> np.ndarray:
    """Percent restoration of a group's profile from Hyper toward Normal.

    Ports analyze_group_true_rescue in legacy/RESCUE_NEW4_1_48_targets_HQ.m::

        D_Hyper       = norm(hyper_profile - normal_profile)
        D_Rescue(c)   = norm(rescued_profiles[:, c] - normal_profile)
        RescuePercent = 100 * (D_Hyper - D_Rescue) / D_Hyper

    Distances are Euclidean over the group's nodes, so a group is scored as a
    profile rather than node-by-node. 100% means the perturbation reaches the
    Normal profile, 0% means no improvement on untreated Hyper, and a negative
    score means it moved the group further away.

    Args:
        hyper_profile: (n_group_nodes,) mean Hyper steady state.
        normal_profile: (n_group_nodes,) mean Normal steady state.
        rescued_profiles: (n_group_nodes, n_strategies) mean perturbed states,
            i.e. the FINAL(rowIdx, :) block of the MATLAB script.

    Returns:
        (n_strategies,) rescue percentages.

    Raises:
        ValueError: if Hyper and Normal are indistinguishable for this group,
            which would make the normalisation meaningless. MATLAB warns and
            skips the group; we refuse explicitly.
    """
    hyper_profile = np.asarray(hyper_profile, dtype=float).ravel()
    normal_profile = np.asarray(normal_profile, dtype=float).ravel()
    rescued_profiles = np.asarray(rescued_profiles, dtype=float)

    if rescued_profiles.ndim != 2:
        raise ValueError(
            f"rescued_profiles must be 2-D (nodes, strategies), got shape "
            f"{rescued_profiles.shape}"
        )
    if not (hyper_profile.shape == normal_profile.shape == (rescued_profiles.shape[0],)):
        raise ValueError(
            "profile shapes disagree: hyper "
            f"{hyper_profile.shape}, normal {normal_profile.shape}, "
            f"rescued {rescued_profiles.shape}"
        )

    d_hyper = float(np.linalg.norm(hyper_profile - normal_profile))
    if d_hyper < np.finfo(float).eps:
        raise ValueError(
            "Hyper and Normal group profiles are indistinguishable; "
            "rescue percentage is undefined for this group"
        )

    d_rescue = np.linalg.norm(
        rescued_profiles - normal_profile[:, None], axis=0
    )
    return 100.0 * (d_hyper - d_rescue) / d_hyper


# N_TOP_TRUE_RESCUE in RESCUE_NEW4_1_48_targets_HQ.m.
DEFAULT_TOP_N = 20


@dataclass(frozen=True)
class GroupRescueRanking:
    """Top-N rescue strategies for one biological group, best first."""

    group: str
    strategy_labels: list[str]     # length min(top_n, n_strategies)
    rescue_percent: np.ndarray     # same length, descending
    distance_to_normal: np.ndarray  # D_Rescue for each listed strategy
    used_nodes: list[str]
    missing_nodes: list[str]
    distance_hyper_to_normal: float


def rank_group_rescue(
    group: str,
    group_nodes: list[str],
    node_names: list[str],
    hyper_mean: np.ndarray,
    normal_mean: np.ndarray,
    final_states: np.ndarray,
    strategy_labels: list[str],
    top_n: int = DEFAULT_TOP_N,
) -> GroupRescueRanking:
    """Rank perturbations by TRUE rescue for one biological group.

    Ports the ranking half of analyze_group_true_rescue. Group membership is
    resolved by exact name match, as MATLAB's `ismember` does — nodes that are
    not in the network are reported in `missing_nodes` rather than silently
    dropped.

    Args:
        group: display name for the group.
        group_nodes: node names making up the group.
        node_names: the network's node names, indexing the arrays below.
        hyper_mean: (n_nodes,) mean Hyper steady state.
        normal_mean: (n_nodes,) mean Normal steady state.
        final_states: (n_nodes, n_strategies), the MATLAB FINAL matrix.
        strategy_labels: one label per column of `final_states`.
        top_n: how many strategies to keep.
    """
    if final_states.shape[1] != len(strategy_labels):
        raise ValueError(
            f"{final_states.shape[1]} strategy columns but "
            f"{len(strategy_labels)} labels"
        )

    index_of = {name: i for i, name in enumerate(node_names)}
    used_nodes = [n for n in group_nodes if n in index_of]
    missing_nodes = [n for n in group_nodes if n not in index_of]
    if not used_nodes:
        raise ValueError(f"group {group!r} has no nodes present in the network")

    rows = [index_of[n] for n in used_nodes]
    hyper_profile = np.asarray(hyper_mean, dtype=float)[rows]
    normal_profile = np.asarray(normal_mean, dtype=float)[rows]
    rescued_profiles = np.asarray(final_states, dtype=float)[rows, :]

    pct = true_rescue_percent(hyper_profile, normal_profile, rescued_profiles)
    d_rescue = np.linalg.norm(rescued_profiles - normal_profile[:, None], axis=0)
    d_hyper = float(np.linalg.norm(hyper_profile - normal_profile))

    keep = min(top_n, pct.size)
    order = np.argsort(-pct, kind="stable")[:keep]

    return GroupRescueRanking(
        group=group,
        strategy_labels=[strategy_labels[i] for i in order],
        rescue_percent=pct[order],
        distance_to_normal=d_rescue[order],
        used_nodes=used_nodes,
        missing_nodes=missing_nodes,
        distance_hyper_to_normal=d_hyper,
    )
