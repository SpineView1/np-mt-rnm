"""Topology metrics (manuscript Figure 2; drawn as TOPO_STATS.png by run_paper_figures.py)."""
from pathlib import Path

import pandas as pd

from np_mt_rnm.network import load_network
from np_mt_rnm.topology import compute_topology

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    net = load_network(ROOT / "data" / "MT_PRIMARY4_1.xlsx")
    metrics = compute_topology(net)
    rows = [
        {
            "node": name,
            "signed_out_degree": metrics.signed_out_degree[name],
            "out_activating": metrics.out_activating[name],
            "out_inhibiting": metrics.out_inhibiting[name],
            "betweenness": metrics.betweenness[name],
            "harmonic_closeness": metrics.harmonic_closeness[name],
        }
        for name in metrics.node_names
    ]
    pd.DataFrame(rows).to_csv(ROOT / "results" / "tables" / "topology_metrics.csv", index=False)


if __name__ == "__main__":
    main()
