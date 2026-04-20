"""Functional-category groupings for the 147 NP-MT network nodes.

Authoritative source: legacy/NP_MT_RNM_FALSIFY4_1.m lines 19–45.
The MATLAB `group_categories` cell array is transcribed verbatim below,
preserving the exact node names from data/MT_PRIMARY4_1.xlsx.

Some nodes appear in multiple categories (e.g., HO-1, GPX1, SOD1, SOD2, CAT
appear in both "metabolic" and "oxidative_proteostasis"). This is intentional
— categories are display groupings, not a strict taxonomy. Use
`get_primary_category(node)` for the single canonical category (first-listed).

The ``"other"`` category collects nodes present in the Excel but not placed in
any MATLAB ``group_categories`` panel. This preserves MATLAB's silence — these
nodes exist in the network but weren't given a display group by the biologist.
No biological assignments have been invented.
"""
from __future__ import annotations

# Short machine-readable IDs → paper display labels (MATLAB's exact strings).
CATEGORY_LABELS: dict[str, str] = {
    "cytokines_chemokines_proteases": "Cytokines, chemokines, proteases & others",
    "growth_factor":                  "Growth factors",
    "transcription_factor":           "Transcription Factors",
    "ion_channel":                    "Ion channels & related",
    "metabolic":                      "Metabolic & related",
    "ecm_matrix":                     "ECM anabolism & phenotype markers",
    "mechanosensor":                  "Mechanical stimuli & its receptors",
    "cell_fate":                      "Cell survival, apoptosis & mitophagy/DNA-damage",
    "oxidative_proteostasis":         "Oxidative-stress defense & proteostasis",
    "mapk":                           "MAPK & stress-activated kinases",
    "rho_cytoskeletal":               "Rho GTPases, cytoskeletal & Hippo regulators",
    "other":                          "Other / unclassified",
}

# Canonical order for figure panel layout (matches MATLAB iteration order).
CATEGORY_ORDER: tuple[str, ...] = (
    "mechanosensor",
    "ion_channel",
    "rho_cytoskeletal",
    "mapk",
    "metabolic",
    "oxidative_proteostasis",
    "transcription_factor",
    "growth_factor",
    "ecm_matrix",
    "cytokines_chemokines_proteases",
    "cell_fate",
    "other",
)

# Transcribed verbatim from MATLAB group_categories. Preserve exact node names.
_CATEGORY_TO_NODES: dict[str, list[str]] = {
    "cytokines_chemokines_proteases": [
        "TNF", "IL6", "IL1β", "IL8", "CCL2", "CXCL1", "CXCL3",
        "ADAMTS4/5", "MMP1", "MMP3", "MMP13",
    ],
    "growth_factor": [
        "TGFβ", "VEGF", "IGF1", "BMP2", "CCN2", "GDF5", "FGF2", "FGF18",
        "Wnt3a", "Wnt5a",
    ],
    "transcription_factor": [
        "CREB", "HIF-1α", "HIF-2α", "NF-κB", "AP-1", "FOXO", "SOX9",
        "NFAT", "RUNX2", "YAP/TAZ", "MRTF-A", "NRF2", "HSF1", "TonEBP",
        "ELK1", "PPARγ", "CITED2",
    ],
    "ion_channel": [
        "Ca2+os", "Ca2+su", "CaMKII", "PKC-E", "PKC-M", "PLCγ-M", "PLCγ-E",
        "CaN", "IP3", "PLA2", "AQP1", "AQP5",
    ],
    "metabolic": [
        "LKB1", "NAD+", "AMPK", "mTORC1", "mTORC2", "SIRT1",
        "PI3K-M", "PI3K-E", "PIP3-M", "PIP3-E", "PDK1-M", "PDK1-E",
        "AKT1-M", "AKT1-E", "GSK3B", "ULK1", "PTEN", "PLD2", "PGE2", "COX-2",
        "CAT", "GPX1", "SOD1", "SOD2", "HO-1", "PHD2", "VHL", "Rheb",
        "NutD", "MitD",
    ],
    "ecm_matrix": [
        "COL2A1", "COL1A1", "COL10A1", "ACAN", "TIMP3",
    ],
    "mechanosensor": [
        "Hypo", "NL", "HL",
        "α5β1-FN", "α5β1-Fs", "αVβ3", "αVβ6",
        "SDC4-E", "SDC4-M", "TRPV4", "PIEZO1",
    ],
    "cell_fate": [
        "Bcl2", "BAX", "CASP3", "CASP9", "BNIP3", "GADD45", "DRP1", "MOMP",
    ],
    "oxidative_proteostasis": [
        "HO-1", "GPX1", "SOD1", "SOD2", "CAT", "HSP70", "HSP27", "ROS",
    ],
    "mapk": [
        "RAS-M", "RAS-E", "RAF-M", "RAF-E", "MEK-M", "MEK-E",
        "ERK-M", "ERK-E", "MKK3/6", "MKK4/7", "JNK", "p38", "P38", "RSK", "TAK1",
    ],
    "rho_cytoskeletal": [
        "RhoA-M", "RhoA-E", "RAC1-M", "RAC1-E", "CDC42",
        "ROCK-M", "ROCK-E", "PAK1", "PKN1",
        "FAK-M", "FAK-E", "MST1/2", "LATS1/2",
    ],
    "other": [
        "ANXA1", "BAD", "DAG", "IKK", "iNOS", "IκBα",
        "LIMK1", "MAP3K", "SMIT1", "SOX5", "SOX6", "TAUT",
        "β-catenin",
    ],
}


# node_name → list of category_ids it belongs to, in CATEGORY_ORDER order.
NODE_CATEGORIES: dict[str, list[str]] = {}
for cat_id in CATEGORY_ORDER:
    for node in _CATEGORY_TO_NODES[cat_id]:
        NODE_CATEGORIES.setdefault(node, []).append(cat_id)


def nodes_in_category(category_id: str) -> list[str]:
    """Return the list of nodes belonging to a category (in paper's order)."""
    if category_id not in _CATEGORY_TO_NODES:
        raise KeyError(f"unknown category_id: {category_id!r}")
    return list(_CATEGORY_TO_NODES[category_id])


def get_primary_category(node_name: str) -> str:
    """Return the single canonical category for a node (first-listed per CATEGORY_ORDER).

    Used for 1-to-1 node coloring in the webapp network view.
    """
    cats = NODE_CATEGORIES.get(node_name)
    if not cats:
        raise KeyError(f"node {node_name!r} has no category assignment")
    return cats[0]


def categories_for(node_name: str) -> list[str]:
    """Return all categories a node belongs to, ordered per CATEGORY_ORDER."""
    return list(NODE_CATEGORIES.get(node_name, ()))
