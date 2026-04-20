# Legacy MATLAB source

These MATLAB files are the original implementation used to produce the paper's
numerical results (Workineh & Noailly, 2026). They are preserved unmodified
as the authoritative reference for the Python port in `../np_mt_rnm/`.

## Running the original MATLAB code

Requires MATLAB R2024a or later with the Parallel Computing Toolbox and
Statistics Toolbox.

```matlab
cd legacy/
NP_MT_RNM_FALSIFY4_1         % baseline simulations + falsification (Figs 4–5, 9)
NP_MT_RNM_FSA4_1             % regime-structured transition paths (Figs 6–8, Supp S3–S6)
RESCUE_NEW4_1                % single + dual-node rescue screen (Figs 10–12, Supp S7–S9)
```

## File index

| File | Responsibility |
|---|---|
| `CreateMatrices_new.m` | Loads `../data/MT_PRIMARY4_1.xlsx` → `Mact`, `Minh` adjacency matrices |
| `ODESysFun.m` | SQUADS ODE RHS (reads Mact/Minh from base workspace; not parallel-safe) |
| `ODESysFunS.m` | SQUADS ODE RHS (explicit parameters; parallel-safe; used in `parfor` loops) |
| `NP_MT_RNM_FALSIFY4_1.m` | Baseline simulations (3 regimes × 100 replicates) + falsification benchmark |
| `NP_MT_RNM_FSA4_1.m` | Regime-structured FSA transition paths + 1D sensitivity scans |
| `RESCUE_NEW4_1.m` | Rescue perturbation screen (single + dual-node clamping) |
