# NP-MT-RNM: Mechanotransduction Regulatory Network Model for Nucleus Pulposus Cells

**A systems-level view of nucleus pulposus mechanotransduction and its regulatory control**

Workineh Z. G.<sup>1</sup>, Chemorion F. K.<sup>1</sup>, Noailly J.<sup>1</sup>

<sup>1</sup>Barcelona Centre for New Medical Technologies (BCN MedTech), Department of Engineering, Universitat Pompeu Fabra, Barcelona, 08018, Spain

---

## Abstract

Intervertebral disc degeneration (IDD), a major contributor to chronic low back pain, is strongly influenced by mechanical loading-dependent regulation of nucleus pulposus (NP) cell phenotype. Previous studies have characterized individual mechanosensors and specific signalling pathways involved in NP mechanotransduction. However, the regulatory interactions through which these components integrate mechanical cues into coordinated cellular responses remain incompletely understood. Here, we present an NP-specific, signed, directed regulatory network model (RNM) comprising 146 molecular nodes and 353 literature-supported interactions. The network connects mechanosensors with signalling cascades, metabolic and redox regulators, transcription factors, extracellular matrix (ECM) effectors, inflammatory mediators, and cell-fate pathways. Using the SQUADS semi-quantitative formalism, which accommodates limited kinetic parameter availability, we simulated steady-state responses to hypo-loading (Hypo), normal loading (NL), and hyper-loading (HL). We also analyzed network topology and screened 48 single- and dual-node perturbations under sustained HL. NL maintained high activity of the ECM-associated regulators ACAN, COL2A1, SOX9, and TIMP3, supported by AKT1, FAK, and antioxidant regulation. HL promoted inflammatory and catabolic regulation through mediators such as NF-κB, AP-1, MMP13, and ADAMTS4/5, accompanied by increased ROS, BAX, and CASP3 activity. Hypo reduced anabolic activity while retaining calcium, ROS, AMPK, and apoptotic signalling, suggesting a possible adhesion-dependent stress mechanism. Model falsification against 39 literature-derived molecular directionality benchmarks yielded 37 concordant predictions (94.9%), providing indirect, direction-level support for the RNM. HSF1 and CCN2 showed contradictory responses. Simulated regulatory recovery was module-specific. SOX9 activation achieved 82.8% ECM-profile recovery, increasing to 85.9% with ROS suppression. SOX9 activation combined with NF-κB suppression achieved 57.0% inflammatory/proteolytic-module recovery. FOXO activation with ROS suppression achieved 68.9% oxidative-stress/proteostasis-module recovery. NRF2 activation with PIEZO1 suppression achieved 37.5% cell-survival-module recovery. These predictions identify candidate molecular targets and intervention combinations for anabolic reinforcement, inflammatory suppression, redox protection, and mechanosensory modulation. The RNM provides a mechanistic framework for computational prioritization of candidate therapeutic targets and intervention combinations. Its predictions can guide subsequent preclinical validation in NP-cell and disc-explant models under controlled mechanical loading. Experimental validation, incorporation of time-dependent loading, and calibration of regulatory parameters are necessary to establish the therapeutic relevance of these predictions and improve the model's translational applicability.

**Keywords:** Nucleus pulposus, Mechanotransduction, Regulatory network modeling, Mechanical loading, Intervertebral disc degeneration, Systems biology, Network control

---

## Overview

This repository contains the computational implementation of a literature-curated **regulatory network model (RNM)** of human **nucleus pulposus (NP) cells** in the intervertebral disc (IVD), specifically extended to capture **mechanotransduction** under chronic mechanical loading. The model integrates upstream mechanosensors, cytoskeletal effectors, MAPK / PI3K / Wnt signaling, transcription factors, and downstream ECM / cytokine / oxidative / cell-fate phenotypes into a single dynamical system driven by three exclusive mechanical inputs: **hypo-loading**, **normal loading**, and **hyper-loading**.

The original MATLAB implementation that produced the paper's numerical results is preserved unmodified under `legacy/` (`NP_MT_FALSIFY4_1_39.m`: baselines and 39-benchmark falsification; `RESCUE_NEW4_1_48_targets_HQ.m`: 48-intervention rescue screen). This Python port regenerates every computational figure of the manuscript from a single command, and additionally exposes the model as an **SBML Level 3 Version 2** file (`model.xml`) so any standards-compliant simulator (libRoadRunner, COPASI, tellurium) can drive it.

### Reproducibility against the manuscript

The port is checked against the manuscript's own numbers, not only against biological plausibility:

- **Rescue screen** — every rescue percentage printed in the manuscript's Results (48 interventions, eight modules), the top intervention of every module, the node-resolved ACAN / COL2A1 / TIMP3 values and the whole Supplementary phenotype-panel table are reproduced **to the printed decimal** (`tests/test_paper_reproduction.py`, `tests/test_rescue_top_perturbations.py`).
- **Falsification** — **37/39** benchmarks concordant (94.9%: 6/6 "+", 31/33 "−"), **HSF1** and **CCN2** contradictory, none indeterminate, as reported (`tests/test_falsification_pass_rate.py`).
- **Topology** — same out-degree, betweenness and harmonic-closeness rankings (NF-κB: 24 outgoing edges, 20 activating / 4 inhibiting; mTORC1 highest betweenness; ROS highest harmonic closeness).

This holds because the port draws its random initial states with MATLAB's own generator: numpy's `RandomState(seed)` is the same Mersenne Twister as MATLAB's `rng(seed, 'twister')`, so with the per-replicate seeds used by `RESCUE_NEW4_1_48_targets_HQ.m` (replicate *r* = 1…100 seeded 1001+*r* / 2001+*r* / 3001+*r* for Hypo / Normal / Hyper) the 100 initial states are bit-identical to MATLAB's.

The model is **multistable**, so ensemble means depend on which initial states are drawn: under Normal loading NRF2, SIRT1 and AMPK settle at 0 in about half of random starts and at 1 in the other half, and under Hyper loading about a quarter of starts remain in the anabolic attractor (ACAN, COL2A1 high). `NP_MT_FALSIFY4_1_39.m` and `NP_MT_RNM_FSA4_1.m` draw `rand` inside `parfor` without per-replicate seeds, so the exact bar heights of the baseline and falsification figures, and the Δ values in the manuscript's Supplementary falsification table, vary between MATLAB runs (e.g. the manuscript reports HL means of 0.78 for BAX, CASP3, CASP9, GSK3B and ROS, the seeded port gives 0.74). The port seeds those ensembles deterministically, reproducing the same charts and conclusions rather than one particular unseeded draw. The falsification outcome does not depend on the draw: across 30 independent 100-replicate Normal/Hyper ensemble pairs it was 37/39, with HSF1 and CCN2 contradictory, in all 30; it is also unchanged with the script's `NL = 0.10` Hyper input (see `legacy/README.md`).

### Network at a glance

| Property | Value |
|---|---|
| Nodes (proteins, ions, ECM components, mechanical inputs) | 146 |
| Directed interactions | 353 |
| &nbsp;&nbsp;&nbsp;&nbsp;Activation edges | 277 |
| &nbsp;&nbsp;&nbsp;&nbsp;Inhibition edges | 76 |
| Mechanical loading inputs (boundary species) | 3 (Hypo, NL, HL) |
| Functional categories (visualization groups) | 12 |
| Falsification benchmarks | 39 (6 NL-higher "+", 33 HL-higher "−"), from 51 literature evidence entries |
| Rescue screen interventions | 48 (6 suppression + 6 activation + 36 dual) |

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

2. **Baseline ensembles**: for each regime, 100 replicate ODE solves from random initial conditions $x_0 \sim \mathcal{U}(0, 1)^N$ (MATLAB-identical seeds, see above), integrated over $t \in [0, 100]$ with RK45 (`rtol = 1e-8`, `atol = 1e-10`, `max_step = 0.5`, as `ode45` in the MATLAB code). The steady state is $x(t = 100)$; every replicate reaches $|dx/dt| < 10^{-8}$. Each node's mean carries a 95 % percentile bootstrap CI (10,000 resamples), and the three regimes are compared per node by Kruskal–Wallis and pairwise Wilcoxon rank-sum tests, BH-FDR adjusted across nodes within each test family (the fixed loading inputs are not tested). These p/q values describe variation among randomised initialisations, not biological replicates.

3. **Falsification benchmark**: 39 literature-derived markers, 6 with higher expression/activity under NL-like loading ("+") and 33 under HL-like loading ("−"). For each node $\Delta_i = \bar{x}_i^{\text{NL}} - \bar{x}_i^{\text{HL}}$ with a 95 % percentile bootstrap CI (10,000 resamples). Zero-based criteria, no tolerance band: a "+" benchmark is *concordant* if the CI lower bound is above 0, a "−" benchmark if the upper bound is below 0; the opposite sign is *contradictory*, and a CI that includes 0 is *indeterminate*.

4. **Hyper → Normal rescue screen**: each of the 100 Hyper steady-state replicates is continued with the intervention clamped (SOX9, PPARγ, HIF-1α, NRF2, SIRT1, FOXO → 1; RhoA-E, PIEZO1, ROS, P38, mTORC1, NF-κB → 0), giving 12 single and 36 cross-class dual interventions (48). Module-level rescue is $R_{p,\mathcal{C}} = 100\,(1 - D_{p,\mathcal{C}} / D_{\text{Hyper},\mathcal{C}})$, the percentage of the Hyper-to-Normal Euclidean distance of the module's ensemble-mean profile removed by intervention $p$; the top 20 per module are reported for eight modules. Module node lists are the rescue script's own readout groups (`np_mt_rnm.rescue.RESCUE_GROUPS`). Node-resolved responses use the paired per-replicate $\Delta_{r,i}^p = x^p_{r,i} - x^{\text{Hyper}}_{r,i}$; phenotype-level readouts average an eight-marker anabolic and a ten-marker catabolic panel.

5. **Constrained transition paths** (earlier manuscript versions; not in the current manuscript, kept as an additional analysis): six steps along Hypo → Normal, $(\text{Hypo}, \text{NL}) = (0.35, 0.10) \to (0.10, 0.35)$ with HL = 0.01, and Normal → Hyper, $(\text{NL}, \text{HL}) = (0.35, 0.10) \to (0.10, 0.35)$ with Hypo = 0.01; 100 replicates per step.

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
│   ├── MT_PRIMARY4_1.xlsx          # model network: 146 nodes, 353 interactions
│   ├── NP-MT-RNM_Supplementary_Data_1_Interaction_Corpus.xlsx  # literature corpus behind every interaction
│   └── falsification_benchmark.csv # 39 benchmarks: direction, evidence IDs, references
│
├── np_mt_rnm/                      # installable Python package
│   ├── __init__.py
│   ├── network.py                  # Excel → adjacency matrices
│   ├── ode.py                      # SQUADS ODE RHS
│   ├── simulation.py               # regime presets, replicate runner, clamps
│   ├── transitions.py              # constrained regime-transition paths + run_trajectory
│   ├── falsification.py            # 39-benchmark registry, zero-based bootstrap-CI criteria
│   ├── rescue.py                   # 48-intervention screen, TRUE-rescue metric
│   ├── statistics.py               # bootstrap CIs, Kruskal-Wallis / rank-sum, permutation, BH-FDR
│   ├── topology.py                 # signed out-degree, betweenness, harmonic closeness
│   ├── categories.py               # functional-category map for all 146 nodes
│   ├── paper_figures.py            # the manuscript's figures, one function per figure
│   ├── paper.py                    # manuscript title, abstract, figure captions (verbatim)
│   ├── figures.py                  # transition-path quick-look figures
│   └── sbml_export.py              # libsbml writer (rate rules)
│
├── scripts/                        # reproducibility entry points
│   ├── run_all.py                  # one-shot full reproduction
│   ├── run_baseline.py             # baseline ensembles, bootstrap CIs, nonparametric tests
│   ├── run_topology.py             # out-degree, betweenness, harmonic closeness
│   ├── run_transitions.py          # constrained transition paths (not in current manuscript)
│   ├── run_falsification.py        # 39-benchmark falsification
│   ├── run_rescue.py               # 48-intervention rescue screen
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
| Fig. 2 — topological metrics | `TOPO_STATS.png` |
| Fig. 3 — baseline: ECM, growth factors, TFs, cytokines (95% bootstrap CIs) | `BSL_ECM_GF.png` |
| Fig. 4 — baseline: metabolic, ion channels, oxidative stress, cell fate | `BSL_MET_OXI.png` |
| Fig. 5A–B — falsification: NL/HL profiles and directional sign grid | `NP_MT_FALS_RECT.png` |
| Fig. 5C — falsification forest plot | `NP_MT_FALS.png` |
| Fig. 6A–D — Top-20 rescue: ECM, growth factors, TFs, cytokines | `Top20_ECM_TrueRescue.png`, `Top20_GF_TrueRescue.png`, `Top20_TF_TrueRescue.png`, `Top20_CYT_TrueRescue.png` |
| Fig. 7A–D — Top-20 rescue: metabolic, ion channels, oxidative stress, cell survival | `Top20_MET_TrueRescue.png`, `Top20_ION_TrueRescue.png`, `Top20_OXI_TrueRescue.png`, `Top20_CS_TrueRescue.png` |
| Figs. S2–S3 — node-resolved responses of the two top interventions per module | `NodeResolved_ECM_GF_TF_CYT.png`, `NodeResolved_MET_ION_SURV_OXI.png` |
| (rescue script output) Δ heatmaps with BH-FDR stars, all / single / dual interventions | `heatmaps/heatmap_{all,singles,duals}_<module>.png` |

Figures 1, 8 and 9 and Figure S1 (Cytoscape network drawing, signalling schematic, falsification workflow, curation pipeline) are illustrations, not model output.

The manuscript captions of Figs. 3–5A state ±1 SD error bars; the figures produced by `NP_MT_FALSIFY4_1_39.m` (and the manuscript PNGs) show 95 % bootstrap CIs of the mean. The port draws CIs, as the MATLAB code does, and writes both CIs and SDs to `results/tables/baseline_means_95ci.csv`.

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
# (REGIME_SEEDS gives the MATLAB seeds of RESCUE_NEW4_1_48_targets_HQ.m)
ens = run_replicates(net, REGIME_PRESETS["Normal"], n_reps=100, seed=REGIME_SEEDS["Normal"])
print(f"SOX9 mean = {ens.mean()[net.node_names.index('SOX9')]:.3f}")

# Single rescue intervention: Hyper + clamp(SOX9=1, ROS=0)
result = run_perturbation(
    net, anabolic_up="SOX9", catabolic_down="ROS",
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

- **Species count:** 146 (one per network node).
- **Floating species:** 143 (carry rate rules).
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

39 rows, one per marker, transcribed from the manuscript's benchmark table (51 evidence entries B01–B51):

```
node, class, expected_direction, literature_direction, benchmark_ids, experimental_evidence, reference_tag
ACAN, anabolic, Normal>Hyper, +, B01;B02, NL-like up; HL-like down, zhang2019mechanical
```

`class` is the MATLAB script's internal label: `anabolic` = higher under NL-like loading ("+"), `catabolic` = higher under HL-like loading ("−"); it does not mean every marker is biologically anabolic or catabolic. `reference_tag` holds the manuscript's BibTeX keys. Order follows `ANABOLIC_LIST` / `CATABOLIC_LIST` of `NP_MT_FALSIFY4_1_39.m`. Note: the manuscript table uses ID B05 for both TIMP3 rows.

---

## Tests

```bash
pytest
```

Tests encode the paper's biological claims as invariants:

- **Paper reproduction** — every rescue percentage printed in the manuscript, to the printed decimal; top-ranked intervention per module; phenotype-panel table.
- **Falsification** — exactly 37/39 concordant, HSF1 and CCN2 contradictory, none indeterminate; zero-based CI classification.
- **MATLAB RNG** — `matlab_rand` returns MT19937 `init_genrand` draws, as `rng(seed, 'twister'); rand`.
- **Network loader** — 146 nodes, 353 interactions (277 activation, 76 inhibition); binary adjacency; no activator/inhibitor overlap.
- **SQUADS ODE** — hand-computed reference values to 12-digit precision.
- **Regime polarities** (Figs. 3–4) — `ACAN[Normal] > 0.8`; `MMP13[Hyper] > 0.8`; `TNF / IL6` near zero under Normal.
- **Rescue design** — 6 suppression + 6 activation targets, 12 singles + 36 duals in MATLAB's loop order; best intervention per module as in the Supplementary table.
- **Baseline statistics** — BH-FDR, bootstrap mean CIs and Kruskal–Wallis / rank-sum against hand computations.
- **Transition paths** — anabolic markers rise monotonically along Hypo → Normal.
- **SBML export** — `model.xml` validates clean; reloading via tellurium reproduces the Normal regime's SOX9-high / ROS-low signature.

---

## Network version

`data/MT_PRIMARY4_1.xlsx` is the network simulated in the manuscript: 146 nodes and 353 signed interactions, exactly the interactions of the curated literature corpus (`data/NP-MT-RNM_Supplementary_Data_1_Interaction_Corpus.xlsx`, Supplementary Data 1 of the paper). The corpus gives, for every interaction, the supporting references (with DOI/PMID), the experimental context, an evidence review of 37 interactions and a Crossref retraction screen of all cited publications.

An earlier version of the network had 357 interactions: four further interactions of a nutrient-deficiency node (`Hypo→NutD`, `NutD→HIF-1α`, `NutD→MitD`, `NutD→ROS`) that are not part of the curated corpus. They and the `NutD` node were removed. With MATLAB-seeded ensembles this changes five printed rescue values by 0.1–0.2 percentage points and left the falsification result of the earlier 45-rule benchmark (43/45) unchanged.

The SQUADS activation term is implemented with the negative exponent, $e^{-h(\omega - 0.5)}$, exactly as in `legacy/ODESysFunS.m`.

---

## Citation

If you use this code or model, please cite both the paper and the software:

> Workineh Z. G., Chemorion F. K. & Noailly J. (2026). *A systems-level view of nucleus pulposus mechanotransduction and its regulatory control.*

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
