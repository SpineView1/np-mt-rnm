# NP-MT-RNM: Mechanotransduction Regulatory Network Model for Nucleus Pulposus Cells

**A systems-level network model reveals mechanical regulation of nucleus pulposus cell states**

Workineh Z. G.<sup>1</sup>, Chemorion F. K.<sup>1</sup>, Noailly J.<sup>1</sup>

<sup>1</sup>BCN MedTech, Universitat Pompeu Fabra, Barcelona, Spain

---

## Overview

This repository contains the computational implementation of a literature-curated **regulatory network model (RNM)** of human **nucleus pulposus (NP) cells** in the intervertebral disc (IVD), specifically extended to capture **mechanotransduction** under chronic mechanical loading. The model integrates upstream mechanosensors, cytoskeletal effectors, MAPK / PI3K / Wnt signaling, transcription factors, and downstream ECM / cytokine / oxidative / cell-fate phenotypes into a single dynamical system driven by three exclusive mechanical inputs: **hypo-loading**, **normal loading**, and **hyper-loading**.

The original MATLAB implementation that produced the paper's numerical results is preserved unmodified under `legacy/`. This Python port regenerates every computational figure of the manuscript from a single command, and additionally exposes the model as an **SBML Level 3 Version 2** file (`model.xml`) so any standards-compliant simulator (libRoadRunner, COPASI, tellurium) can drive it.

### Reproducibility against the manuscript

The port is checked against the manuscript's own numbers, not only against biological plausibility:

- **Rescue screen** — every percentage printed in the manuscript text and on the bars of the six Top-20 rescue figures is reproduced **to the printed decimal** (`tests/test_paper_reproduction.py`).
- **Falsification** — **43/45** rules satisfied, failing exactly **TonEBP** and **HSF1**, as reported (`tests/test_falsification_pass_rate.py`).
- **Topology** — same out-degree, betweenness and harmonic-closeness rankings (NF-κB: 24 outgoing edges, 20 activating / 4 inhibiting; mTORC1 highest betweenness; ROS highest harmonic closeness).

This holds because the port draws its random initial states with MATLAB's own generator: numpy's `RandomState(seed)` is the same Mersenne Twister as MATLAB's `rng(seed, 'twister')`, so with the per-replicate seeds used by `RESCUE_NEW4_1_final.m` (replicate *r* = 1…100 seeded 1001+*r* / 2001+*r* / 3001+*r* for Hypo / Normal / Hyper) the 100 initial states are bit-identical to MATLAB's.

The model is **multistable**, so ensemble means depend on which initial states are drawn: under Normal loading NRF2, SIRT1 and AMPK settle at 0 in about half of random starts and at 1 in the other half, and under Hyper loading about a quarter of starts remain in the anabolic attractor (ACAN, COL2A1 high). `NP_MT_RNM_FALSIFY4_1.m` and `NP_MT_RNM_FSA4_1.m` draw `rand` inside `parfor` without per-replicate seeds, so the exact bar heights of the baseline, falsification and transition figures vary between MATLAB runs; the port seeds those ensembles deterministically, reproducing the same charts and conclusions (pass/fail outcome, regime patterns) rather than one particular unseeded draw. Across 30 independent 100-replicate ensembles the falsification outcome was 43/45 in 28 and 41/45 (NRF2 and SIRT1 also failing) in 2.

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
| Hypo  | 0.20 | 0.01 | 0.01 |
| Normal | 0.01 | 0.80 | 0.01 |
| Hyper  | 0.01 | 0.01 | 0.80 |

2. **Baseline ensembles**: for each regime, 100 replicate ODE solves from random initial conditions $x_0 \sim \mathcal{U}(0, 1)^N$ (MATLAB-identical seeds, see above), integrated over $t \in [0, 100]$ with RK45 (`rtol = 1e-8`, `atol = 1e-10`, `max_step = 0.5`, as `ode45` in the MATLAB code). The steady state is $x(t = 100)$; every replicate reaches $|dx/dt| < 10^{-8}$.

3. **Falsification benchmark**: 45 rules (17 anabolic, expected Normal > Hyper; 28 catabolic, expected Hyper > Normal). For each node $\Delta = \bar{x}_{\text{Normal}} - \bar{x}_{\text{Hyper}}$ with a 95 % nonparametric bootstrap CI (10,000 resamples). An anabolic rule passes if the CI lower bound exceeds $F_{\text{TOL}} = 0.02$; a catabolic rule if the upper bound is below $-0.02$.

4. **Hyper → Normal rescue screen**: each of the 100 Hyper steady-state replicates is continued with the intervention clamped (SOX9, PPARγ, HIF-1α, NRF2, IκBα → 1; RhoA-E, PIEZO1, PI3K-E, FAK-E, ROS → 0), giving 10 single and 25 cross-class dual interventions. Module-level rescue is $R_{p,\mathcal{C}} = 100\,(1 - D_{p,\mathcal{C}} / D_{\text{Hyper},\mathcal{C}})$, the percentage of the Hyper-to-Normal Euclidean distance of the module's ensemble-mean profile removed by intervention $p$; the top 20 per module are reported. Node-resolved responses use the paired per-replicate $\Delta_{r,i}^p = x^p_{r,i} - x^{\text{Hyper}}_{r,i}$.

5. **Constrained transition paths**: six steps along Hypo → Normal, $(\text{Hypo}, \text{NL}) = (0.35, 0.10) \to (0.10, 0.35)$ with HL = 0.01, and Normal → Hyper, $(\text{NL}, \text{HL}) = (0.35, 0.10) \to (0.10, 0.35)$ with Hypo = 0.01; 100 replicates per step.

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
│   ├── MT_PRIMARY4_1.xlsx          # 147 nodes, 357 edges (network used for the manuscript)
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
│   ├── run_baseline.py             # baseline ensembles (Hypo / Normal / Hyper)
│   ├── run_topology.py             # out-degree, betweenness, harmonic closeness
│   ├── run_transitions.py          # constrained transition paths
│   ├── run_falsification.py        # 45-rule benchmark + pass rate
│   ├── run_rescue.py               # 35-intervention rescue screen
│   ├── run_paper_figures.py        # every manuscript figure → results/figures/paper/
│   ├── build_web_bundle.py         # JSON artifacts for the webapp
│   ├── build_paper_bundle.py       # paper_results.json (webapp "Paper results" tab)
│   └── export_sbml.py              # regenerate model.xml
│
├── results/                        # committed build artifacts
│   ├── figures/paper/              # the manuscript's figures, same file names
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

- Python ≥ 3.11
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

- `results/figures/paper/` — every computational figure of the manuscript, under the manuscript's own file names (below)
- `results/tables/` — one CSV per analysis
- `results/replicates/` — raw steady-state ensembles (NPZ)
- `results/web_bundle/` — JSON artifacts consumed by [np-mt-rnm-web](https://github.com/SpineView1/np-mt-rnm-web)

Full reproduction takes about 4 minutes on 10 cores (dominated by the rescue screen).

| Manuscript figure | File in `results/figures/paper/` |
|---|---|
| Fig. 3 — topological metrics | `TOPO_STATS.png` |
| Fig. 4 — baseline: ECM, growth factors, TFs, cytokines | `BSL_ECM_GF.png` |
| Fig. 5 — baseline: metabolic, ion channels, oxidative stress, cell fate | `BSL_MET_OXI.png` |
| Fig. 6 — ECM transition heatmaps | `ECM_transition.png` |
| Fig. 7 — representative transition trajectories | `Representative_transition_paths.png` |
| Fig. 8 — falsification | `NP_MT_FALS.png` |
| Fig. 9 — ECM rescue ranking | `ECM_rescue.png` |
| Supp. S5 — transition heatmaps | `TF_transition.png`, `GF_transition.png`, `CYT_transition.png`, `OXI_transition.png`, `CSF_transition.png` |
| Supp. S6 — rescue rankings | `GF_rescue.png`, `TF_rescue.png`, `CYT_rescue.png`, `OXI_rescue.png`, `CSF_rescue.png` |
| Supp. S6 — node-resolved rescue responses | `ECM_rescue1.png`, `GF_rescue1.png`, `TF_rescue1.png`, `CYT_rescue1.png`, `OXI_rescue1.png`, `CSF_rescue1.png` |

Figures 1–2 and Supplementary S1 (pathway schematic, Cytoscape network drawing, curation pipeline) are illustrations, not model output.

---

## Usage as a Python library

```python
from np_mt_rnm.network import load_network
from np_mt_rnm.simulation import REGIME_PRESETS, REGIME_SEEDS, run_replicates
from np_mt_rnm.rescue import run_perturbation
from np_mt_rnm.transitions import run_trajectory, HYPO_TO_NORMAL_PATH

# Load topology
net = load_network("data/MT_PRIMARY4_1.xlsx")
print(f"{len(net.node_names)} nodes, "
      f"{int(net.mact.sum())} act + {int(net.minh.sum())} inh edges")

# Baseline ensemble under Normal loading
# (REGIME_SEEDS gives the MATLAB seeds of RESCUE_NEW4_1_final.m)
ens = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=100, seed=REGIME_SEEDS["Normal"])
print(f"SOX9 mean = {ens.mean()[net.node_names.index('SOX9')]:.3f}")

# Single rescue perturbation: Hyper + clamp(SOX9=1, RhoA-E=0)
result = run_perturbation(
    net, anabolic_up="SOX9", catabolic_down="RhoA-E",
    n_reps=100,
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

- **Paper reproduction** — every rescue percentage printed in the manuscript, to the printed decimal; top-ranked intervention per module.
- **Falsification** — exactly 43/45, failing TonEBP and HSF1.
- **MATLAB RNG** — `matlab_rand` returns MT19937 `init_genrand` draws, as `rng(seed, 'twister'); rand`.
- **Network loader** — 147 nodes, 357 edges (281 activation, 76 inhibition); binary adjacency; no activator/inhibitor overlap.
- **SQUADS ODE** — hand-computed reference values to 12-digit precision.
- **Regime polarities** (Figs. 4–5) — `ACAN[Normal] > 0.8`; `MMP13[Hyper] > 0.8`; `TNF / IL6` near zero under Normal.
- **Rescue top perturbations** — node-level response rankings (mean |Δ|) per module.
- **Transition paths** (Figs. 6–7) — anabolic markers rise monotonically along Hypo → Normal.
- **SBML export** — `model.xml` validates clean; reloading via tellurium reproduces the Normal regime's SOX9-high / ROS-low signature.

---

## Network version

`data/MT_PRIMARY4_1.xlsx` is the **357-edge** network that produced the manuscript's figures. A later revision of the spreadsheet (2026-09) removed the four `NutD` edges (`Hypo→NutD`, `NutD→HIF-1α`, `NutD→MitD`, `NutD→ROS`), leaving 353; with it the transcription-factor and oxidative-stress rescue percentages move by 0.1–0.3 points and no longer match the printed values. The repository follows the manuscript.

The SQUADS activation term is implemented with the negative exponent, $e^{-h(\omega - 0.5)}$, exactly as in `legacy/ODESysFunS.m`.

---

## Citation

If you use this code or model, please cite both the paper and the software:

> Workineh Z. G., Chemorion F. K. & Noailly J. (2026). *A systems-level network model reveals mechanical regulation of nucleus pulposus cell states.*

Software citation: see `CITATION.cff`.

---

## Authors

| Author | Role | Affiliation |
|---|---|---|
| Zerihun G. Workineh | Lead author, model design, MATLAB implementation | BCN MedTech, Universitat Pompeu Fabra |
| Francis K. Chemorion | Python port, SBML export, software engineering | BCN MedTech, Universitat Pompeu Fabra |
| Jérôme Noailly | Senior author, supervision | BCN MedTech, Universitat Pompeu Fabra |

---

## License

MIT — see [`LICENSE`](LICENSE).
