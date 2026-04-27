"""Generate the SBML representation of the NP-MT regulatory network.

Outputs:
  model.xml — SBML Level 3 Version 2 model. Each non-boundary node carries a
              RateRule encoding the SQUADS/Mendoza ODE; Hypo / NL / HL are
              boundary species, clampable by the loader.

Usage:
    python scripts/export_sbml.py [--output model.xml] [--h 10] [--gamma 1]
"""
from __future__ import annotations

import argparse
from pathlib import Path

from np_mt_rnm.network import load_network
from np_mt_rnm.sbml_export import export_sbml

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "MT_PRIMARY4_1.xlsx"


def main() -> None:
    parser = argparse.ArgumentParser(description="Export NP-MT-RNM as SBML")
    parser.add_argument(
        "--output", "-o", default=str(ROOT / "model.xml"),
        help="Output path (default: <repo-root>/model.xml)",
    )
    parser.add_argument("--h", type=float, default=10.0, help="Sigmoid gain (default 10)")
    parser.add_argument("--gamma", type=float, default=1.0, help="Decay constant (default 1)")
    args = parser.parse_args()

    net = load_network(DATA)
    print(f"Network: {len(net.node_names)} nodes, "
          f"{int(net.mact.sum())} activation edges, "
          f"{int(net.minh.sum())} inhibition edges")

    ok = export_sbml(net, args.output, h=args.h, gamma=args.gamma)
    if not ok:
        raise SystemExit("SBML export failed")


if __name__ == "__main__":
    main()
