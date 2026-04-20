"""One-shot reproduction: runs every analysis script in dependency order.

Use this for a clean rebuild of results/. Each step can also be run
individually via its own script.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STEPS = [
    "scripts/run_baseline.py",
    "scripts/run_topology.py",
    "scripts/run_transitions.py",
    "scripts/run_falsification.py",
    "scripts/run_rescue.py",
    "scripts/build_web_bundle.py",
]


def main() -> None:
    for step in STEPS:
        print(f"\n===== {step} =====", flush=True)
        result = subprocess.run([sys.executable, str(ROOT / step)], cwd=ROOT)
        if result.returncode != 0:
            print(f"FAILED: {step} (exit code {result.returncode})")
            sys.exit(result.returncode)
    print("\nAll steps complete. Results in results/.")


if __name__ == "__main__":
    main()
