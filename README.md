# NP-MT-RNM: Mechanotransduction Regulatory Network Model for Nucleus Pulposus Cells

**A systems-level network model reveals how mechanical loading organizes regulatory states and transitions in nucleus pulposus cells**

Workineh Z. G.<sup>1</sup>, Noailly J.<sup>1</sup>

<sup>1</sup>BCN MedTech, Universitat Pompeu Fabra, Barcelona, Spain

---

## Overview

This repository contains the computational implementation of a literature-curated **regulatory network model (RNM)** of human **nucleus pulposus (NP) cells** in the intervertebral disc (IVD), specifically extended to capture **mechanotransduction** under chronic mechanical loading. The model integrates upstream mechanosensors, cytoskeletal effectors, MAPK / PI3K / Wnt signaling, transcription factors, and downstream ECM / cytokine / oxidative / cell-fate phenotypes into a single dynamical system driven by three exclusive mechanical inputs: **hypo-loading**, **normal loading**, and **hyper-loading**.

The original MATLAB implementation that produced the paper's numerical results is preserved unmodified under `legacy/`. This Python port reproduces every paper figure (Figs 3–12 + Supp S3–S9) from a single command, and additionally exposes the model as an **SBML Level 3 Version 2** file (`model.xml`) so any standards-compliant simulator (libRoadRunner, COPASI, tellurium) can drive it.

### Network at a glance

| Property | Value |
|---|---|
| Nodes (proteins, ions, ECM components, mechanical inputs) | 147 |
| Directed interactions | 357 |
| &nbsp;&nbsp;&nbsp;&nbsp;Activation edges | 281 |
| &nbsp;&nbsp;&nbsp;&nbsp;Inhibition edges | 76 |
| Mechanical loading inputs (boundary species) | 3 (Hypo, NL, HL) |
| Functional categories (visualization groups) | 12 |
| Falsification benchmark rules | 45 (17 anabolic + 28 catabolic) |
| Rescue screen perturbations | 35 (5 KD + 5 UP + 25 dual) |

---

## Mathematical framework

### SQUADS / Mendoza ODE system

The static knowledge-based topology is converted into a **semi-quantitative dynamical system** following the SQUADS formalism — a parameterised case of the framework of [Mendoza & Xenarios (2006)](https://doi.org/10.1186/1742-4682-3-13). Each node activation evolves according to an ODE that integrates upstream activating and inhibiting signals through a sigmoidal transfer function.

#### State equation

For each protein node $x_n$:

$$\frac{dx_n}{dt} = \frac{-e^{0.5 h} + e^{-h(\omega_n - 0.5)}}{(1 - e^{0.5 h})(1 + e^{-h(\omega_n - 0.5)})} - \gamma\, x_n$$

where:

- $h$ — sigmoid **gain**. Larger $h$ pushes the system toward Boolean (ON/OFF) behaviour.
- $\omega_n \in [0, 1]$ — **aggregated regulatory input** combining all upstream activators and inhibitors of node $n$.
- $\gamma$ — linear **decay** rate.

The first term is a sigmoid activation function mapping $\omega_n \in [0, 1]$ to an activation level in $[0, 1]$.

#### Aggregated regulatory input $\omega_n$

The sign of $\omega_n$ depends on which classes of upstream regulators are present. With uniform interaction weights $\alpha = \beta = 1$:

**Case (i) — only activators** present ($\{x_{nk}^a\}$ non-empty, no inhibitors):

$$\omega_n = \frac{k_a + 1}{k_a} \cdot \frac{\sum_k x_{nk}^a}{1 + \sum_k x_{nk}^a}$$

**Case (ii) — only inhibitors** present ($\{x_{nl}^i\}$ non-empty, no activators):

$$\omega_n = 1 - \frac{k_i + 1}{k_i} \cdot \frac{\sum_l x_{nl}^i}{1 + \sum_l x_{nl}^i}$$

**Case (iii) — both activators and inhibitors:**

$$\omega_n = \left(\frac{k_a + 1}{k_a} \cdot \frac{\sum_k x_{nk}^a}{1 + \sum_k x_{nk}^a}\right) \cdot \left(1 - \frac{k_i + 1}{k_i} \cdot \frac{\sum_l x_{nl}^i}{1 + \sum_l x_{nl}^i}\right)$$

where $k_a$ and $k_i$ count the activators and inhibitors of node $n$, respectively. This formulation keeps $\omega_n \in [0, 1]$, preserving the Boolean asymptotic limit at high gain.

**Case (iv) — pure inputs** (no regulators, e.g. Hypo, NL, HL): $\omega_n$ is undefined; the node is held at its clamp value as a boundary species.

#### Default parameter values

| Parameter | Symbol | Value | Rationale |
|---|---|---|---|
| Sigmoid gain | $h$ | 10 | Intermediate steepness; matches paper |
| Decay constant | $\gamma$ | 1 | Uniform across all nodes |
| Activation weights | $\alpha$ | 1 | Uniform (no per-edge sensitivity data) |
| Inhibition weights | $\beta$ | 1 | Uniform |

### Simulation protocol

1. **Three mechanical loading regimes** define the boundary inputs:

| Regime | Hypo | NL | HL |
|---|---|---|---|
| Hypo  | 0.80 | 0.01 | 0.01 |
| Normal | 0.01 | 0.80 | 0.01 |
| Hyper  | 0.01 | 0.01 | 0.80 |

2. **Baseline ensembles** (Figs 4–5): for each regime, 100 replicate ODE solves from independent random initial conditions $x_0 \sim \mathcal{U}(0, 1)^N$, integrated over $t \in [0, 100]$. The steady state is taken as $x(t = 100)$.

3. **Falsification benchmark** (Fig 9): 45 hand-curated rules of the form *"node X is expected to go up/down between Normal and Hyper"*, scored against the model with paired bootstrap CIs (10,000 resamples, BH-FDR for multiple testing).

4. **Rescue screen** (Figs 10–12 + Supp S7–S9): starting from the Hyper steady-state ensemble, 35 perturbations evaluate single (5 catabolic ↓ + 5 anabolic ↑) and dual ($5 \times 5$) clamps; ranked by per-category mean $|\Delta|$ relative to the unperturbed Hyper baseline.

5. **Regime transitions** (Figs 6–8 + Supp S3–S6): constrained gradient paths along Hypo→Normal and Normal→Hyper, computed as a sequence of intermediate boundary inputs with replicate-paired solves at each step.

---

## Repository structure

```
np-mt-rnm/
├── README.md                       # this file
├── model.xml                       # SBML Level 3 V2 export (auto-generated)
├── pyproject.toml
├── requirements.txt
├── CITATION.cff
│
├── data/
│   ├── MT_PRIMARY4_1.xlsx          # 147-node × 4-col edge list
│   └── falsification_benchmark.csv # 45 curated rules
│
├── np_mt_rnm/                      # installable Python package
│   ├── __init__.py
│   ├── network.py                  # Excel → adjacency matrices
│   ├── ode.py                      # SQUADS ODE RHS
│   ├── simulation.py               # regime presets, replicate runner, clamps
│   ├── transitions.py              # constrained regime-transition paths + run_trajectory
│   ├── falsification.py            # 45-rule benchmark + bootstrap CI
│   ├── rescue.py                   # 35-perturbation screen
│   ├── statistics.py               # Cohen's d, bootstrap, permutation, BH-FDR
│   ├── topology.py                 # signed out-degree, betweenness, harmonic closeness
│   ├── categories.py               # functional-category map for all 147 nodes
│   ├── figures.py                  # matplotlib code matching paper's visual style
│   └── sbml_export.py              # libsbml writer (rate rules)
│
├── scripts/                        # reproducibility entry points
│   ├── run_all.py                  # one-shot full reproduction
│   ├── run_baseline.py             # Figs 3, 4, 5
│   ├── run_topology.py             # Fig 3
│   ├── run_transitions.py          # Figs 6, 7, 8 + Supp S3–S6
│   ├── run_falsification.py        # Fig 9 + pass rate
│   ├── run_rescue.py               # Figs 10, 11, 12 + Supp S7–S9
│   ├── build_web_bundle.py         # JSON artifacts for the webapp
│   └── export_sbml.py              # regenerate model.xml
│
├── results/                        # committed build artifacts
│   ├── figures/                    # Figs 3–12 + Supp S3–S9 (PNG)
│   ├── tables/                     # one CSV per analysis
│   ├── replicates/                 # raw steady-state ensembles (NPZ)
│   └── web_bundle/                 # JSON consumed by np-mt-rnm-web
│
├── legacy/                         # original MATLAB code, preserved unmodified
└── tests/                          # biological-invariance tests
```

---

## Installation

### Requirements

- Python ≥ 3.10
- NumPy, SciPy, pandas, openpyxl, matplotlib
- joblib (parallel ensemble runs)
- python-libsbml (SBML export)
- tellurium (optional; only needed to simulate the SBML through libRoadRunner)

### Setup

```bash
git clone https://github.com/SpineView1/np-mt-rnm.git
cd np-mt-rnm
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

### Reproduce all paper figures

```bash
python scripts/run_all.py
```

Outputs land in:

- `results/figures/` — Figures 3–12 + Supp S3–S9 (PNG)
- `results/tables/` — one CSV per analysis
- `results/replicates/` — raw steady-state ensembles (NPZ)
- `results/web_bundle/` — JSON artifacts consumed by [np-mt-rnm-web](https://github.com/SpineView1/np-mt-rnm-web)

Full reproduction takes 5–15 minutes on 8 cores (dominated by the rescue screen at ~2.5 min + transitions at ~30 s).

---

## Usage as a Python library

```python
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, run_replicates
from np_mt_rnm.rescue import run_perturbation
from np_mt_rnm.transitions import run_trajectory, HYPO_TO_NORMAL_PATH

# Load topology
net = load_network("data/MT_PRIMARY4_1.xlsx")
print(f"{len(net.node_names)} nodes, "
      f"{int(net.mact.sum())} act + {int(net.minh.sum())} inh edges")

# Baseline ensemble under Normal loading
ens = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=100, seed=42)
print(f"SOX9 mean = {ens.mean()[net.node_names.index('SOX9')]:.3f}")

# Single rescue perturbation: Hyper + clamp(SOX9=1, RhoA-E=0)
result = run_perturbation(
    net, anabolic_up="SOX9", catabolic_down="RhoA-E",
    n_reps=100, seed=42,
)
print(f"|Δ| over all nodes = {abs(result.mean_delta).mean():.3f}")

# Continuous trajectory across a Normal→Hyper switch at t=10
traj = run_trajectory(
    net, REGIME_PRESETS["Normal"], REGIME_PRESETS["Hyper"],
    t_switch=10.0, t_end=30.0, n_t=31, n_reps=10, seed=42, n_jobs=1,
)
print(f"trajectory shape (n_t, n_nodes) = {traj.mean.shape}")
```

---

## SBML model

The repository ships a pre-generated `model.xml` (SBML Level 3 Version 2) encoding the complete SQUADS / Mendoza ODE system as **rate rules** rather than per-reaction kinetics:

```text
d[species]/dt = sigmoid(omega(activators, inhibitors); h) - gamma * [species]
```

- **Species count:** 147 (one per network node).
- **Floating species:** 144 (carry rate rules).
- **Boundary species:** 3 — `Hypo`, `NL`, `HL`. Marked `boundaryCondition=true` so any SBML loader can clamp them via `runner["Hypo"] = value` (or equivalent), reproducing the three loading regimes without modifying the model.
- **Global parameters:** `h = 10`, `gamma = 1`.
- **Validation:** validates clean (`document.checkConsistency()` returns 0 errors).

### Regenerate the SBML

```bash
python scripts/export_sbml.py --output model.xml
```

### Simulate the SBML directly (e.g. via tellurium)

```python
import tellurium as te

r = te.loadSBMLModel("model.xml")
r["Hypo"] = 0.01
r["NL"]   = 0.80
r["HL"]   = 0.01
result = r.simulate(0, 100, 51)
sox9_final = result["[SOX9]"][-1]
print(f"SOX9 final under Normal regime: {sox9_final:.3f}")
```

The accompanying webapp [np-mt-rnm-web](https://github.com/SpineView1/np-mt-rnm-web) consumes this same SBML.

---

## Data format

### Network edge list (`data/MT_PRIMARY4_1.xlsx`)

One row per node (target):

| Nodes | Activators | Inhibitors | Stimuli |
|---|---|---|---|
| ACAN | SOX9, SMAD2/3 | NF-κB | ACAN |
| Hypo | NOTHING | NOTHING | Hypo |

`Activators` and `Inhibitors` are comma-separated upstream node names. `NOTHING` (or empty) indicates the node has no regulators of that type. The `Stimuli` column duplicates the node name and is used by the legacy MATLAB clamp-builder; in the Python port, the three pure-input nodes (Hypo, NL, HL) are the only species that act as boundary conditions.

### Falsification benchmark (`data/falsification_benchmark.csv`)

45 rules. Each row:

```
node, class, expected, reference_tag
```

`class` ∈ {anabolic, catabolic}; `expected` ∈ {up, down}; `reference_tag` is a short citation key linking back to the supporting literature (full bibliography in the paper).

---

## Tests

```bash
pytest
```

Tests encode the paper's biological claims as invariants:

- **Network loader** — 147 nodes, 357 edges; binary adjacency; no activator/inhibitor overlap.
- **SQUADS ODE** — hand-computed reference values to 12-digit precision.
- **Regime polarities** (Figs 4–5) — `ACAN[Normal] > 0.8`; `MMP13[Hyper] > 0.8`; `TNF / IL6` near zero under Normal.
- **Falsification pass rate** (Fig 9) — ≥ 88 % concordance (paper: 95.6 %; Python port: 91.1 % due to RNG differences vs. MATLAB; 2 borderline-CI failures).
- **Rescue top perturbations** (Figs 10B, 12B) — top ECM rescuer involves SOX9; top TF rescuer involves NRF2.
- **Transition paths** (Fig 8) — anabolic markers rise monotonically along Hypo → Normal.
- **SBML export** — `model.xml` validates clean; reloading via tellurium reproduces the Normal regime's SOX9-high / ROS-low signature.

---

## Known paper-text corrections

The Python port surfaced two typos / inconsistencies in the submitted paper text:

1. **Edge count** — paper Section 2.1 + abstract say *356 edges*; the Excel actually contains **357**. Section 3.1 (which already says 357) is correct.
2. **SQUADS formula sign** — paper Section 2.2 prints `exp(h·(ω − 0.5))` in the activation term; the MATLAB source (and the physically correct form, which is what produced the paper's numerical results) uses `exp(−h·(ω − 0.5))`.

Both should be corrected in the manuscript before submission.

---

## Citation

If you use this code or model, please cite both the paper and the software:

> Workineh Z. G. & Noailly J. (2026). *A systems-level network model reveals how mechanical loading organizes regulatory states and transitions in nucleus pulposus cells.*

Software citation: see `CITATION.cff`.

---

## Authors

| Author | Role | Affiliation |
|---|---|---|
| Zerihun G. Workineh | Lead author, model design, MATLAB implementation | BCN MedTech, Universitat Pompeu Fabra |
| Jérôme Noailly | Senior author, supervision | BCN MedTech, Universitat Pompeu Fabra |
| Francis K. Chemorion | Python port, SBML export, software engineering | BCN MedTech, Universitat Pompeu Fabra |

---

## License

MIT — see [`LICENSE`](LICENSE).
