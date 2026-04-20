# np-mt-rnm

Python implementation of the NP cell mechanotransduction regulatory network
model accompanying Workineh & Noailly (2026), *A systems-level network model
reveals how mechanical loading organizes regulatory states and transitions
in nucleus pulposus cells*.

## Quick reproduction

```bash
git clone https://github.com/SpineView1/np-mt-rnm.git
cd np-mt-rnm
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run_all.py
```

Outputs land in `results/figures/`, `results/tables/`, and
`results/web_bundle/`.

## Repo layout

- `np_mt_rnm/` — installable Python package
- `data/` — network edge list + falsification benchmark
- `scripts/` — reproducibility entry points
- `results/` — committed build artifacts (figures, tables, JSON bundle)
- `legacy/` — original MATLAB code (reference only, preserved unmodified)
- `tests/` — biological-invariance tests encoding the paper's claims

## Citing

See `CITATION.cff`.
