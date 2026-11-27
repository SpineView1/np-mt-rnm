# Legacy MATLAB source

These MATLAB files are the original implementation used to produce the
manuscript's numerical results (Workineh, Chemorion & Noailly, 2026). They are
preserved unmodified as the authoritative reference for the Python port in
`../np_mt_rnm/`.

## Running the original MATLAB code

Requires MATLAB R2024a or later with the Parallel Computing Toolbox and
Statistics and Machine Learning Toolbox.

```matlab
cd legacy/
NP_MT_FALSIFY4_1_39            % baselines with 95% bootstrap CIs, KW/rank-sum + BH-FDR,
                               % 39-benchmark falsification (manuscript Figs 3-5)
RESCUE_NEW4_1_48_targets_HQ    % 6 suppression x 6 activation targets: 12 single + 36 dual
                               % node clamps, TRUE-rescue rankings (Figs 6-7, Supp. Figs S2-S3)
NP_MT_RNM_FSA4_1               % constrained transition paths (not in the current manuscript)
```

## File index

| File | Responsibility |
|---|---|
| `CreateMatrices_new.m` | Loads `../data/MT_PRIMARY4_1.xlsx` → `Mact`, `Minh` adjacency matrices |
| `ODESysFun.m` | SQUADS ODE RHS (reads Mact/Minh from base workspace; not parallel-safe) |
| `ODESysFunS.m` | SQUADS ODE RHS (explicit parameters; parallel-safe; used in `parfor` loops) |
| `NP_MT_FALSIFY4_1_39.m` | Baselines (3 regimes × 100 replicates) with bootstrap mean CIs and nonparametric tests; falsification against 39 benchmarks (6 NL-higher, 33 HL-higher) with zero-based 95% bootstrap CI criteria |
| `RESCUE_NEW4_1_48_targets_HQ.m` | Hyper→Normal rescue screen: suppress RhoA-E, PIEZO1, ROS, P38, mTORC1, NF-κB (→0); activate SOX9, PPARγ, HIF-1α, NRF2, SIRT1, FOXO (→1); 48 strategies |
| `NP_MT_RNM_FSA4_1.m` | Regime-structured FSA transition paths + 1D sensitivity scans |

## Revision history

The previous versions, `NP_MT_RNM_FALSIFY4_1.m` (45 benchmarks, ±0.02 tolerance
band) and `RESCUE_NEW4_1_final.m` (5 × 5 targets, 35 strategies), are in the
git history of this directory; the files were renamed to the new script names.

## Known discrepancy

`NP_MT_FALSIFY4_1_39.m` sets `NL_hype = 0.10` for the Hyper regime, whereas the
manuscript Methods and `RESCUE_NEW4_1_48_targets_HQ.m` use `NL = 0.01`. The
script itself flags this and keeps 0.10; the Python port follows the
manuscript (0.01). With seeded ensembles the falsification outcome (37/39,
HSF1 and CCN2 contradictory) is the same.
