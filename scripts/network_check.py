"""Does a target's position in the model's semantic network predict model-human (dis)agreement?

Exploratory follow-up to the pilot's network analysis. For the configuration that won the
ManyLanguages discovery search (outputs/discovery_grid.csv, row 1):

 1. Encode every stimulus word (targets and both cue sets) alone, at the layer chosen for the
    similarity covariate (best |r| with human priming on discovery, from outputs/layer_profile.csv).
 2. Neighbourhood measures per target, relative to all stimulus words:
      density = mean cosine to its 10 nearest neighbours
      degree  = number of words above the 95th percentile of all pairwise cosines
 3. Test, on each half, whether density/degree add to human priming beyond
    [dS, dSim, log_freq, length, semd], whether they moderate the dS-human link (interaction),
    and whether the dS-human correlation differs across density terciles.

The discovery half is the primary test; the confirmation half checks it. Nothing here selects a
configuration.

Usage: python scripts/network_check.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, get_split, load_human  # noqa: E402
from compute_priming import layer_embeddings, is_masked, pick_device  # noqa: E402

OUT = ROOT / "outputs"
K = 10
COVARS = ["dSim", "log_freq", "length", "semd"]


def z(v):
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std()


def fit(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return beta, 1 - resid.var() / y.var()


def added(y, base, extra):
    """R^2 gain, standardized beta and p for the LAST column of `extra` added to `base`."""
    _, r2b = fit(y, base)
    full = np.column_stack([base, extra])
    beta, r2f = fit(y, full)
    n, k = len(y), full.shape[1] + 1
    f = (r2f - r2b) / ((1 - r2f) / (n - k))
    return r2f - r2b, beta[-1], stats.f.sf(f, 1, n - k)


def main():
    grid = pd.read_csv(OUT / "discovery_grid.csv")
    best = grid.iloc[0]
    model, fmt = best["model"], best["format"]
    print(f"Configuration: {model} | {fmt} (discovery r = {best['r_dS']:.3f})")

    prof = pd.read_csv(OUT / "layer_profile.csv")
    pm = prof[(prof.model == model) & (prof.half == "discovery")]
    layer = int(pm.loc[pm["r_human"].abs().idxmax(), "layer"])
    print(f"Layer (chosen on discovery): {layer}")

    slug = "".join(c if c.isalnum() or c in ".-" else "_" for c in model)
    run = pd.read_csv(OUT / "surprisal" / f"{slug}__{fmt}.csv")
    words = sorted(set(run.target) | set(run.cue_related) | set(run.cue_unrelated))
    vecs = layer_embeddings(model, pick_device(), words, is_masked(model))
    V = np.stack([vecs[w][layer] for w in words])
    V = V / np.linalg.norm(V, axis=1, keepdims=True)
    S = V @ V.T
    np.fill_diagonal(S, -np.inf)
    thr = np.percentile(S[np.isfinite(S)], 95)
    idx = {w: i for i, w in enumerate(words)}
    nb = pd.DataFrame({
        "target": run.target,
        "density": [np.sort(S[idx[t]])[-K:].mean() for t in run.target],
        "degree": [(S[idx[t]] > thr).sum() for t in run.target],
    })
    print(f"{len(words)} stimulus words; degree threshold (95th pct cosine) = {thr:.3f}")

    d = (run[["target", "dS", f"dSim_L{layer:02d}"]].rename(columns={f"dSim_L{layer:02d}": "dSim"})
         .merge(load_human(), on="target").merge(nb, on="target")
         .merge(get_split()[["target", "half"]], on="target").dropna())
    d["log_degree"] = np.log1p(d["degree"])

    print("\nCorrelation of neighbourhood measures with other variables (all items):")
    for m in ["density", "log_degree"]:
        print(f"  {m}: with avgZ_prime r={stats.spearmanr(d[m], d.avgZ_prime)[0]:.3f}; "
              f"with dS r={stats.spearmanr(d[m], d.dS)[0]:.3f}; with log_freq r={stats.spearmanr(d[m], d.log_freq)[0]:.3f}")

    rows = []
    for half in ["discovery", "confirmation"]:
        h = d[d.half == half]
        y = z(h.avgZ_prime)
        base = np.column_stack([z(h[c]) for c in ["dS"] + COVARS])
        for m in ["density", "log_degree"]:
            dR2, b, p = added(y, base, z(h[m]))
            inter = z(h["dS"]) * z(h[m])
            dR2i, bi, pi = added(y, np.column_stack([base, z(h[m])]), inter)
            rows.append(dict(half=half, measure=m, n=len(h), main_dR2=dR2, main_beta=b, main_p=p,
                             interaction_dR2=dR2i, interaction_beta=bi, interaction_p=pi))
    res = pd.DataFrame(rows)
    print("\n=== Added to [dS, dSim, log_freq, length, semd]: main effect and dS x measure interaction ===")
    print(res.round(4).to_string(index=False))
    res.to_csv(OUT / "network_check.csv", index=False)

    d["tercile"] = pd.qcut(d["density"], 3, labels=["low", "mid", "high"])
    print("\n=== Spearman r(dS, human priming) by density tercile ===")
    t = []
    for half in ["discovery", "confirmation"]:
        for lab, g in d[d.half == half].groupby("tercile", observed=True):
            t.append(dict(half=half, density=lab, n=len(g), r=stats.spearmanr(g.dS, g.avgZ_prime)[0]))
    print(pd.DataFrame(t).pivot(index="density", columns="half", values="r").round(3).to_string())

    hubs = d.sort_values("degree", ascending=False).head(10)[["target", "degree", "density", "avgZ_prime"]]
    print("\nTop 10 hub targets:")
    print(hubs.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
