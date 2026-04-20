# np-mt-rnm

Python implementation of the nucleus pulposus (NP) cell mechanotransduction
regulatory network model accompanying:

> Workineh Z. G. & Noailly J. (2026). *A systems-level network model reveals
> how mechanical loading organizes regulatory states and transitions in
> nucleus pulposus cells.*

**147 molecular nodes · 357 experimentally supported interactions · three
mechanical loading regimes (hypo, normal, hyper) · 45-rule falsification
benchmark · 35-perturbation rescue screen.**

The original MATLAB implementation used to generate the paper's numerical
results is preserved unmodified under `legacy/`. This Python port reproduces
every paper figure from a single command.

## Reproduce the paper figures

```bash
git clone https://github.com/SpineView1/np-mt-rnm.git
cd np-mt-rnm
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run_all.py
```

Outputs land in:

- `results/figures/` — Figures 3–12 + Supp S3–S9 (PNG)
- `results/tables/` — one CSV per analysis
- `results/replicates/` — raw steady-state ensembles (NPZ)
- `results/web_bundle/` — JSON artifacts consumed by the webapp
  ([np-mt-rnm-web](https://github.com/SpineView1/np-mt-rnm-web))

Full reproduction takes 5–15 minutes on 8 cores (dominated by the rescue
screen at ~2.5 min + transitions at ~30 s).

## Repo layout

```
np_mt_rnm/           # installable Python package
├── network.py       # edge-list loader (ports CreateMatrices_new.m)
├── ode.py           # SQUADS ODE RHS (ports ODESysFunS.m)
├── simulation.py    # regime presets, replicate runner, clamps
├── transitions.py   # constrained regime transition paths (ports FSA script)
├── falsification.py # 45-rule benchmark + bootstrap CI
├── rescue.py        # 35-perturbation screen (ports RESCUE script)
├── statistics.py    # Cohen's d, bootstrap, permutation, BH-FDR
├── topology.py      # signed out-degree, betweenness, harmonic closeness
├── categories.py    # functional-category map for all 147 nodes
└── figures.py       # matplotlib code matching paper's visual style

data/
├── MT_PRIMARY4_1.xlsx         # network edge list (147 nodes × 4 cols)
└── falsification_benchmark.csv # 45 curated rules (17 anabolic + 28 catabolic)

scripts/             # reproducibility entry points
├── run_all.py       # one-shot full reproduction
├── run_baseline.py  # Figs 3, 4, 5
├── run_topology.py  # Fig 3
├── run_transitions.py # Figs 6, 7, 8 + Supp S3–S6
├── run_falsification.py # Fig 9 + pass rate
├── run_rescue.py    # Figs 10, 11, 12 + Supp S7–S9
└── build_web_bundle.py  # JSON artifacts for webapp

results/             # committed build artifacts (figures, tables, web_bundle)
legacy/              # original MATLAB code, preserved unmodified
tests/               # biological-invariance tests encoding the paper's claims
docs/                # design documentation
```

## Tests

```bash
pytest
```

47 tests encode the paper's biological claims as invariants:

- **Network loader** — 147 nodes, 357 edges; binary adjacency; no activator/inhibitor overlap.
- **SQUADS ODE** — hand-computed reference values to 12-digit precision.
- **Regime polarities** (Figs 4–5) — ACAN Normal > 0.8; MMP13 Hyper > 0.8; TNF/IL6 Normal near zero.
- **Falsification pass rate** (Fig 9) — ≥ 88% concordance (paper: 95.6%; Python port: 91.1% due to RNG differences; see commit history for details).
- **Rescue top perturbations** (Figs 10B, 12B) — top ECM rescuer involves SOX9; top TF rescuer involves NRF2.
- **Transition paths** (Fig 8) — anabolic markers rise monotonically along Hypo→Normal.

## Known paper-text corrections

The Python port surfaced two typos/inconsistencies in the submitted paper text:

1. **Edge count**: paper Section 2.1 says 356 edges; Excel actually contains
   **357**. Section 3.1 (which already says 357) is correct.
2. **SQUADS formula sign**: paper Section 2.2 displays `exp(h·(ω−0.5))` in the
   activation term; the MATLAB source (and physically-correct form) uses
   `exp(−h·(ω−0.5))`. The numerical results in the paper used the MATLAB formula.

Both should be corrected before submission.

## License

MIT — see `LICENSE`.

## Citation

See `CITATION.cff`.
