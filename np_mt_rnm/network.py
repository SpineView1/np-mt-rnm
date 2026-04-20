"""Network loader. Ports legacy/CreateMatrices_new.m.

The MATLAB function reads an Excel file with columns
[Nodes, Activators, Inhibitors, Stimuli] and builds signed adjacency
matrices Mact (147x147) and Minh (147x147). This module reproduces that
behavior using pandas + numpy.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_NULL_TOKENS = {"", "NOTHING", "NONE", "NA", "N/A", "NULL", "-"}


@dataclass(frozen=True)
class Network:
    """Loaded NP-MT regulatory network."""

    node_names: list[str]
    mact: np.ndarray  # shape (n, n), binary
    minh: np.ndarray  # shape (n, n), binary
    stimuli_names: list[str]


def _tokenize(cell: object) -> list[str]:
    """Split a cell like 'Hypo, TRPV4' or 'NOTHING' or '' into trimmed tokens.

    Empty strings and null placeholders are dropped. Case is preserved but
    comparisons upstream should be case-insensitive.
    """
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return []
    tokens = [tok.strip() for tok in str(cell).split(",")]
    return [t for t in tokens if t and t.upper() not in _NULL_TOKENS]


def load_network(xlsx_path: Path | str) -> Network:
    """Load the edge list into a `Network`.

    Ports `CreateMatrices_new.m`. Column semantics:
      Nodes      — node name (row = target)
      Activators — comma-separated node names that activate this target
      Inhibitors — comma-separated node names that inhibit this target
      Stimuli    — the stimulus label for this row (blank for non-input nodes)
    """
    df = pd.read_excel(xlsx_path, sheet_name=0)
    required_cols = {"Nodes", "Activators", "Inhibitors", "Stimuli"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Excel missing required columns: {sorted(missing)}")

    node_names = [str(n).strip() for n in df["Nodes"].tolist()]
    n = len(node_names)
    name_to_idx = {name.lower(): i for i, name in enumerate(node_names)}

    stimuli_names = [
        str(s).strip()
        for s in df["Stimuli"].tolist()
        if pd.notna(s) and str(s).strip().upper() not in _NULL_TOKENS
    ]

    mact = np.zeros((n, n), dtype=np.int8)
    minh = np.zeros((n, n), dtype=np.int8)

    for i, (act_cell, inh_cell) in enumerate(zip(df["Activators"], df["Inhibitors"])):
        for token in _tokenize(act_cell):
            j = name_to_idx.get(token.lower())
            if j is None:
                raise ValueError(
                    f"Row {i} ({node_names[i]}): activator '{token}' not in node list"
                )
            mact[i, j] = 1
        for token in _tokenize(inh_cell):
            j = name_to_idx.get(token.lower())
            if j is None:
                raise ValueError(
                    f"Row {i} ({node_names[i]}): inhibitor '{token}' not in node list"
                )
            minh[i, j] = 1

    return Network(
        node_names=node_names,
        mact=mact,
        minh=minh,
        stimuli_names=stimuli_names,
    )
