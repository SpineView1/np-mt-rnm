# NP-MT-RNM Paper Artifact — Design Spec

**Date:** 2026-04-20
**Author:** Francis Chemorion (UPF), with Claude
**Paper:** Workineh & Noailly, *A systems-level network model reveals how mechanical loading organizes regulatory states and transitions in nucleus pulposus cells* (2026)
**Submission target:** ~2026-04-27

## 1. Context

The biologist (Zerihun Workineh) has MATLAB code implementing a regulatory network model (RNM) of NP cell mechanotransduction: 147 nodes, 356 signed edges, three loading regimes (Hypo, Normal, Hyper), a falsification benchmark (45 rules, 95% pass rate), constrained transition-path analysis, and a rescue perturbation screen.

The paper is being submitted in ~1 week and requires:

1. A **public, citable GitHub repository** containing the model code — the paper artifact.
2. A **Django webapp** for interactive exploration, deployed to a UPF server via Docker — for lab and collaborator use.

**Existing MATLAB code** (`/Users/kiptengwer/Documents/ZERIHUN/`):
- `CreateMatrices_new.m` — loads `MT_PRIMARY4_1.xlsx` edge list → `Mact`, `Minh` adjacency matrices.
- `ODESysFun.m` / `ODESysFunS.m` — SQUADS ODE right-hand side (serial / parallel-safe variants).
- `NP_MT_RNM_FALSIFY4_1.m` — baseline simulations (3 regimes × 100 replicates) + falsification.
- `NP_MT_RNM_FSA4_1.m` — regime-structured FSA transition paths.
- `RESCUE_NEW4_1.m` — single + dual-node rescue perturbation screen.

**Reference repos** (prior work, structural templates only):
- `github.com/kchemorion/RNMUpdate` — earlier Python port (has clean `rnm/` package layout).
- `github.com/SpineView1/RNM` — earlier Django webapp (has tab-partial template layout).

## 2. Decisions

| Decision | Value |
|---|---|
| Repo strategy | Two fresh repos, no extending old ones |
| Model repo | `SpineView1/np-mt-rnm` (public, MIT) |
| Webapp repo | `SpineView1/np-mt-rnm-web` (public, MIT) |
| Roles | Model repo = paper artifact (Zenodo-minted after v1); webapp = lab-internal interactive tool |
| Port strategy | Full Python port; original MATLAB preserved in `legacy/` |
| Python version | 3.11 |
| Webapp framework | Django + Cytoscape.js + Chart.js |
| Webapp compute | Live 1-replicate ODE for Simulate + custom Rescue; heavy workflows (100 reps, bootstrap, full rescue screen) precomputed and served as static JSON |
| Deployment | Docker image for UPF internal servers (no cloud CI required) |
| Webapp tabs | Overview · Network · Simulate · Falsification · Transitions · Rescue · Downloads (7) |
| Simulate tab layout | Top-strip controls, split network + bars below |
| Rescue tab layout | Side-by-side: published heatmap \| design-your-own |
| Figure fidelity | Matplotlib code visually matches paper figures (colors, fonts, layout) |
| Testing | Biological-invariance tests (regime polarities, falsification pass rate, rescue top perturbations) + numerical stability (`max\|dx/dt\| < 1e-8`). **No cached MATLAB reference fixtures** — Python scripts run end-to-end from scratch, same as biologist runs MATLAB. |

## 3. Model Repo Architecture

```
SpineView1/np-mt-rnm/
├── np_mt_rnm/                           # installable Python package
│   ├── __init__.py
│   ├── network.py                       # load edge list → Mact, Minh, NodeNames
│   ├── ode.py                           # SQUADS RHS + clamped wrapper + solver config
│   ├── simulation.py                    # baseline/replicate runner, regime helpers
│   ├── rescue.py                        # single/dual perturbation enumeration + screen
│   ├── transitions.py                   # FSA regime-structured paths + 1D scans
│   ├── falsification.py                 # bootstrap CI, FTOL rules, benchmark loader
│   ├── statistics.py                    # Cohen's d, permutation tests, BH-FDR, bootstrap
│   ├── categories.py                    # functional-category mapping (single source of truth)
│   └── figures.py                       # matplotlib code tuned to paper visual style
├── data/
│   ├── MT_PRIMARY4_1.xlsx               # edge list (147 nodes × 4 cols)
│   └── falsification_benchmark.csv      # 45 curated rules from paper Table S2
├── scripts/
│   ├── run_baseline.py                  # Figs 4–5 + replicate data
│   ├── run_falsification.py             # Fig 9 + benchmark pass rate
│   ├── run_transitions.py               # Figs 6–8 + supp S3–S6
│   ├── run_rescue.py                    # Figs 10–12 + supp S7–S9
│   ├── build_web_bundle.py              # CSVs → web_bundle/*.json
│   └── run_all.py                       # full reproduction one-shot
├── results/                             # committed build artifacts
│   ├── figures/                         # PNGs matching paper figures
│   ├── tables/                          # CSVs
│   ├── replicates/                      # .npz raw steady-state ensembles
│   └── web_bundle/                      # JSON for webapp consumption
├── legacy/                              # original MATLAB, untouched
│   ├── NP_MT_RNM_FALSIFY4_1.m
│   ├── NP_MT_RNM_FSA4_1.m
│   ├── RESCUE_NEW4_1.m
│   ├── ODESysFun.m, ODESysFunS.m, CreateMatrices_new.m
│   └── README.md                        # how to run the original MATLAB
├── tests/
│   ├── test_network.py                  # edge-list loader correctness
│   ├── test_ode.py                      # SQUADS math at fixed inputs
│   ├── test_regime_polarities.py        # Normal/Hypo/Hyper produce expected polarities
│   ├── test_falsification_pass_rate.py  # full falsification ≥ 95% pass rate
│   ├── test_rescue_top_perturbations.py # top ECM rescuer = SOX9; top TF rescuer = NRF2
│   └── test_ode_steady_state.py         # max|dx/dt| < 1e-8 at t=100
├── pyproject.toml
├── requirements.txt
├── CITATION.cff
├── LICENSE
└── README.md
```

## 4. Webapp Repo Architecture

```
SpineView1/np-mt-rnm-web/
├── np_mt_rnm_web/                       # Django project
│   ├── settings.py, urls.py, wsgi.py
├── rnm_app/                             # Django app
│   ├── views.py                         # tab views + /api/simulate, /api/rescue/custom
│   ├── urls.py
│   ├── templates/rnm/
│   │   ├── base.html                    # 7-tab nav (reference SpineView1/RNM style)
│   │   ├── partials/overview.html
│   │   ├── partials/network.html        # Cytoscape.js container
│   │   ├── partials/simulate.html       # Layout B: top strip + split network/bars
│   │   ├── partials/falsification.html
│   │   ├── partials/transitions.html
│   │   ├── partials/rescue.html         # Layout B: published | design-your-own
│   │   └── partials/downloads.html
│   └── static/rnm/
│       ├── js/cytoscape.min.js, chart.min.js
│       ├── js/simulate.js, rescue.js, network.js
│       └── css/rnm.css
├── precomputed/                         # synced from model repo's results/web_bundle/
├── Dockerfile                           # gunicorn + whitenoise
├── docker-compose.yml
├── requirements.txt                     # pins np-mt-rnm @ git+...@v0.1.0
├── LICENSE
└── README.md
```

**Model/webapp coupling:** webapp pins the model package via `np-mt-rnm @ git+https://github.com/SpineView1/np-mt-rnm.git@v0.1.0`. Updating the model = bump tag + rebuild Docker image. Webapp never touches the MATLAB code or raw Excel; it imports the Python package and reads `web_bundle/*.json`.

## 5. Data Flow

### One-time reproducibility run (model repo)

```
MT_PRIMARY4_1.xlsx  →  network.load_network()  →  (Mact, Minh, NodeNames)
                                                         │
                                                         ▼
                        ode.squads_rhs(t, X, Mact, Minh, clamps, h=10, γ=1)
                                                         │
                                                         ▼
                           simulation.run_replicates(regime, n_reps=100)
                                                         │
                          ┌──────────────────┼──────────────────┐
                          ▼                  ▼                  ▼
                falsification.screen  transitions.fsa_scan  rescue.screen
                          │                  │                  │
                          └──────────────────┼──────────────────┘
                                             ▼
                                figures.py + tables (matplotlib, pandas)
                                             ▼
                  results/figures/*.png + results/tables/*.csv + results/web_bundle/*.json
```

### Runtime (webapp live endpoint)

```
Browser (Simulate)                      Django view
──────────────────                      ───────────
[user drags slider]                     POST /api/simulate
[JS debounces 150ms]  ──HTTP JSON──►   {regime: {...}, clamps: {"SOX9": 1}}
                                                │
                                                ▼
                                      simulation.run_single_replicate(...)
                                                │  (~0.6–1.0 s)
                                                ▼
                                      {node_activations: {...147}}
Chart.js + Cytoscape.js  ◄──JSON──    renders bars + recolors network
```

Stateless POST per simulation. No database, no sessions, no Celery.

## 6. Numerical Methods

- **ODE formalism:** SQUADS (Sánchez-Rodríguez style), exactly as in `ODESysFun.m`. Per-node activation `X_i ∈ [0,1]`. `h = 10`, `γ = 1`, `α = β = 1`. Combined activation/inhibition via `ω = A·(1−I)` when both present.
- **Solver:** `scipy.integrate.solve_ivp(method='RK45', rtol=1e-8, atol=1e-10, max_step=0.5)`.
- **Time window:** `t ∈ [0, 100]`.
- **Steady-state criterion:** `max|dx/dt| < 1e-8` at `t=100` (assertion in tests).
- **Clamped nodes:** Hypo/NL/HL inputs (and user-supplied clamps) override the ODE. Matches MATLAB wrapper `ODESysFun_Clamped_safe3` exactly — inside the RHS, clamped indices have their derivative forced to `0` (so they can't drift), and an assertion after `solve_ivp` confirms the final `X` values at clamped indices equal the clamp values to within `atol`.
- **Regime presets:** Hypo = (0.20, 0.01, 0.01); Normal = (0.01, 0.80, 0.01); Hyper = (0.01, 0.01, 0.80).
- **Replicates:** random `X_0 ∈ [0,1]^147` per replicate; 100 replicates per regime for paper-matching results.
- **Parallelization:** `joblib.Parallel(n_jobs=-1)` for replicate loops and rescue screen. Not shared-memory — each worker re-initializes `Mact`, `Minh`.

## 7. API Contract

### `POST /api/simulate`

Request:
```json
{
  "regime": {"Hypo": 0.01, "NL": 0.80, "HL": 0.01},
  "clamps": {"SOX9": 1.0, "NF-kB": 0.0}
}
```

Response (~0.6–1.0 s):
```json
{
  "node_activations": {"SOX9": 1.0, "ACAN": 0.87, "...": 147},
  "converged": true,
  "max_abs_derivative": 2.3e-9,
  "t_final": 100.0
}
```

### `POST /api/rescue/custom`

Request:
```json
{
  "anabolic_up": "SOX9",
  "catabolic_down": "RhoA-E",
  "n_replicates": 20
}
```

Server caps `n_replicates` via the Django setting `RNM_RESCUE_LIVE_MAX_REPLICATES` (default `20`, ~15–25 s response). UI displays "quick preview — paper uses 100 replicates".

Response:
```json
{
  "perturbation_label": "SOX9↑ + RhoA-E↓",
  "mean_delta_per_node": {"...": "147 entries grouped by category"},
  "balance_score": 0.42,
  "distance_to_normal": 0.18
}
```

### `GET /api/network`

Cached from `precomputed/network.json` (synced at Docker build time from the model repo's `results/web_bundle/network.json`). Returns nodes + edges + functional-category labels + layout coordinates.

## 8. Precomputed Artifacts

```
results/
├── figures/                             # PNGs — visual layout matches paper (colors, fonts, panel arrangement)
│   ├── fig2_network.png
│   ├── fig3_topology_metrics.png
│   ├── fig4_baseline_ecm_gf_tf_cyto.png
│   ├── fig5_baseline_metab_ion_ox_apop.png
│   ├── fig6_transition_ecm.png
│   ├── fig7_transition_tf.png
│   ├── fig8_transition_trajectories.png
│   ├── fig9_falsification_forest.png
│   ├── fig10_rescue_ecm.png
│   ├── fig11_rescue_growth_factors.png
│   ├── fig12_rescue_transcription.png
│   └── supp/S3.png ... S9.png
├── tables/
│   ├── baseline_summary.csv
│   ├── falsification_rules.csv
│   ├── rescue_summary.csv
│   └── rescue_per_node_deltas.csv
├── replicates/
│   ├── baseline_hypo.npz
│   ├── baseline_normal.npz
│   └── baseline_hyper.npz
└── web_bundle/
    ├── network.json
    ├── baseline.json
    ├── falsification.json
    ├── rescue.json
    └── transitions.json
```

## 9. Testing Strategy

The biologist's workflow is to run the MATLAB scripts from scratch and validate figures by eye. The Python port follows the same workflow. No cached MATLAB reference fixtures are used.

**Automated tests encode the paper's own claims as invariants:**

- `test_network.py` — edge-list loader produces 147 nodes; edge count matches the Excel (reconcile 356 vs 357 discrepancy between paper Sections 2.1 and 3.1).
- `test_ode.py` — SQUADS RHS at a fixed input vector gives the expected numerical output: two canned test cases (one single-activator node, one combined activator+inhibitor node) computed from the SQUADS formula by hand in a Jupyter notebook committed alongside the test, confirming sign and magnitude of `dx/dt`.
- `test_ode_steady_state.py` — all regimes converge with `max|dx/dt| < 1e-8` at `t=100`.
- `test_regime_polarities.py` — per Figs 4–5: Normal `ACAN > 0.8` and `COL10A1 ≈ 0`; Hyper `COL10A1 > 0.8` and `ACAN < 0.3`; MMP13 near zero under Normal, high under Hyper; TNF/IL6 near zero under Normal.
- `test_falsification_pass_rate.py` — full falsification run asserts ≥ 95% pass rate (paper: 43/45 = 95.6%).
- `test_rescue_top_perturbations.py` — top-ranked ECM rescuer involves SOX9 activation (Fig 10B); top-ranked transcription-factor rescuer involves NRF2 (Fig 12B).

**Validation gate before ship:** reviewer opens each regenerated figure next to the paper PDF and confirms visual match. Same process the biologist uses in MATLAB today.

## 10. Build Sequence

1. **Model repo skeleton.** Scaffold `np-mt-rnm` (pyproject, MIT, pre-commit). Port `CreateMatrices_new.m` → `network.py`. Port `ODESysFunS.m` → `ode.py`. Implement `simulation.run_single_replicate(regime, clamps)`. *Done when:* single ODE solve for Normal regime converges with `max|dx/dt| < 1e-8`.
2. **Baseline reproduction (Figs 3/4/5).** Multi-replicate runner (100 reps × 3 regimes, joblib). `categories.py` functional-group map. Matplotlib figure code matching paper style. *Done when:* Figs 3–5 regenerated visually match paper; `test_regime_polarities.py` green.
3. **Transitions / FSA (Figs 6/7/8 + Supp S3–S6).** Port `NP_MT_RNM_FSA4_1.m` → `transitions.py`. *Done when:* transition heatmaps + trajectory plots reproduced.
4. **Falsification (Fig 9).** Port falsification section → `falsification.py`. Bootstrap CI (10k resamples), FTOL rules, forest plot. *Done when:* Fig 9 A+B reproduced; `test_falsification_pass_rate.py` green.
5. **Rescue screen (Figs 10/11/12 + Supp S7–S9).** Port `RESCUE_NEW4_1.m` → `rescue.py`. 35 perturbations × 100 reps. `scripts/run_all.py` produces full `results/` + `web_bundle/` end-to-end. *Done when:* Figs 10–12 reproduced; `test_rescue_top_perturbations.py` green.
6. **Webapp skeleton.** Django project + 7-tab shell with stub partials. Dockerfile (gunicorn + whitenoise). Pip-install `np-mt-rnm` from GitHub tag. *Done when:* `docker compose up` serves the 7-tab page with `web_bundle/` mounted.
7. **Network + Simulate tabs (live).** Network tab renders Cytoscape.js from `web_bundle/network.json`. Simulate tab (Layout B) wired to `POST /api/simulate`. *Done when:* dragging a slider or adding a clamp updates the bar chart + network coloring in ≤1 s.
8. **Falsification + Transitions + Downloads tabs (display-only).** Static renders from `web_bundle/`. *Done when:* forest plot + transition heatmaps match the paper in-browser.
9. **Rescue tab (Layout B).** Left = precomputed heatmap + top-10 ranking from `web_bundle/rescue.json`. Right = "design your own" calling `POST /api/rescue/custom`. *Done when:* published heatmap matches Figs 10A/11A/12A; custom combo returns deltas in < 30 s.
10. **Ship polish.** `CITATION.cff`, README with clone-and-reproduce + screenshot. Tag `v0.1.0` on model repo. Build webapp Docker image for UPF registry. *Done when:* fresh clone + `docker compose up` works end-to-end.

**Dependencies.** Steps 1→2→3→4→5 strictly sequential in the model repo. Step 6 can start after Step 2 (needs network + baseline JSON). Step 7 after Step 2. Step 8 after Step 4. Step 9 after Step 5. Step 10 last.

**Downscope ladder** (drop in this order only if a checkpoint slips):

1. Supplementary heatmaps (S3–S9) — point to paper PDF.
2. Rescue "design your own" live mode → display-only.
3. Network tab → static PNG instead of Cytoscape.js.
4. Simulate tab → display-only with canned example traces.

## 11. Open Items

- **Edge count discrepancy.** Paper Section 2.1 says 356 edges; Section 3.1 says 357. Resolve against `MT_PRIMARY4_1.xlsx` when loader is built. Update paper if needed.
- **Supplementary `NP-MT-EdgeList-Topology` file.** Paper mentions a structured edge-list file with literature annotations per edge. Need to confirm whether this is `MT_PRIMARY4_1.xlsx` or a separate artifact. If separate, that file goes in `data/` too.
- **Parallel ODE solve reproducibility.** Each joblib worker seeds its RNG independently. Document the seeding protocol in `simulation.py` so figures are regenerable bit-for-bit across machines with the same `n_jobs`.
- **UPF Docker registry details.** Hostname, auth method, push procedure — clarify before Step 10.
- **Zenodo release.** Out of scope for v1; list in post-submission checklist.
