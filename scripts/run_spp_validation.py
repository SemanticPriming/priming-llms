"""Run the SPP external validation for configs chosen ONLY from the ManyLanguages discovery grid.

Rule (fixed here, before any SPP result exists): the TOP_N configurations by discovery
Spearman r in outputs/discovery_grid.csv. Run scripts/analyze.py first to write that file.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TOP_N = 8

grid = pd.read_csv(ROOT / "outputs" / "discovery_grid.csv").head(TOP_N)
print("SPP validation configs (top", TOP_N, "on ManyLanguages discovery):")
print(grid[["config", "r_dS"]].round(3).to_string(index=False), flush=True)
for kind in ["first", "other"]:
    for model, g in grid.groupby("model"):
        subprocess.run([sys.executable, str(ROOT / "scripts" / "compute_priming.py"),
                        "--models", model, "--formats", *g["format"],
                        "--items-file", str(ROOT / "data" / "spp" / f"spp_items_{kind}.csv"),
                        "--out-dir", str(ROOT / "outputs" / f"spp_{kind}")], check=True)
