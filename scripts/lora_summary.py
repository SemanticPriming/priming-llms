"""Combine the multi-seed LoRA runs: real vs shuffled-label gain over the untrained model.

Reads outputs/lora/*__summary.csv written by lora_align.py (filter: rank 4, 3 epochs). For each
evaluation set it reports, across seeds, the mean gain (r_lora - r_untrained) for real and shuffled
labels, and a paired test of real minus shuffled (the real gain beyond what capacity alone gives).
The loaded sets: confirmation (primary), and the original SPP outcomes (external).

Usage: python scripts/lora_summary.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
frames = [pd.read_csv(f) for f in (ROOT / "outputs" / "lora").glob("*__summary.csv")]
d = pd.concat(frames, ignore_index=True)
d = d[(d["rank"] == 4) & (d["epochs"] == 3)] if "epochs" in d else d.iloc[0:0]
if d.empty:
    raise SystemExit("No multi-seed LoRA summaries found yet.")
d = d[d["set"].isin(["confirmation", "spp_first", "spp_other"])]
d["variant"] = d["model"].str.split("/").str[-1] + " | " + d["format"]
d["label"] = d["set"] + " | " + d["outcome"]

rows = []
for (var, lab), g in d.groupby(["variant", "label"]):
    real = g[~g.shuffled].set_index("seed")
    shuf = g[g.shuffled].set_index("seed")
    seeds = sorted(set(real.index) & set(shuf.index))
    if len(seeds) < 2:
        continue
    gr, gs = real.loc[seeds, "gain"].to_numpy(), shuf.loc[seeds, "gain"].to_numpy()
    diff = gr - gs
    t = stats.ttest_1samp(diff, 0) if len(seeds) > 2 else None
    rows.append(dict(variant=var, evaluation=lab, seeds=len(seeds), r_untrained=real.loc[seeds, "r_untrained"].mean(),
                     gain_real=gr.mean(), gain_real_sd=gr.std(ddof=1), gain_shuffled=gs.mean(),
                     real_minus_shuffled=diff.mean(),
                     ci_lo=diff.mean() - stats.t.ppf(.975, len(seeds) - 1) * diff.std(ddof=1) / np.sqrt(len(seeds)),
                     ci_hi=diff.mean() + stats.t.ppf(.975, len(seeds) - 1) * diff.std(ddof=1) / np.sqrt(len(seeds)),
                     p=getattr(t, "pvalue", np.nan)))
out = pd.DataFrame(rows)
print(out.round(3).to_string(index=False))
out.to_csv(ROOT / "outputs" / "lora_multiseed_summary.csv", index=False)
