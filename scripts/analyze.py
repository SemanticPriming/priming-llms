"""Discovery/confirmation analysis of model priming (dS) vs human priming (avgZ_prime).

Reads every outputs/surprisal/*.csv written by compute_priming.py, then:

 A. LAYER PROFILE   For each model, correlate each layer's static similarity (dSim_L) with
                    (i) human priming and (ii) the model's own dS. The layer used as the
                    similarity covariate is chosen on the DISCOVERY half only.
 B. DISCOVERY GRID  Every model x format scored by Spearman r(dS, avgZ_prime) on discovery.
                    The winner is picked by SELECTION_RULE, fixed before confirmation is read.
 C. CONFIRMATION    The winner is evaluated once on the confirmation half; we also test whether
                    dS adds to [dSim at the chosen layer, frequency, length, semantic diversity].
 D. FEW-SHOT / FRAME  Descriptive tables: confirmation r by frame and by number of shots.
 E. READOUT CEILING A cross-validated ridge fit on discovery (all dS + similarity + covariates),
                    scored on confirmation. SUPERVISED: report as an upper bound, not as the
                    model's own priming.

Usage: python scripts/analyze.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, SEED, get_split, load_human  # noqa: E402

SELECTION_RULE = "max discovery Spearman r(dS, avgZ_prime)"  # fixed before confirmation is read
OUT = ROOT / "outputs"
COVARS = ["log_freq", "length", "semd"]


def spearman_ci(x, y, n_boot=1000, seed=SEED):
    rng = np.random.default_rng(seed)
    x, y = np.asarray(x), np.asarray(y)
    r = stats.spearmanr(x, y)[0]
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        boots.append(stats.spearmanr(x[i], y[i])[0])
    return r, *np.percentile(boots, [2.5, 97.5])


def z(v):
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std()


def ols_r2(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta, 1 - (y - X @ beta).var() / y.var()


def incremental(d, covars):
    y = z(d["avgZ_prime"])
    C = np.column_stack([z(d[c]) for c in covars])
    _, r2_base = ols_r2(y, C)
    beta, r2_full = ols_r2(y, np.column_stack([C, z(d["dS"])]))
    n, k = len(y), C.shape[1] + 2
    f = (r2_full - r2_base) / ((1 - r2_full) / (n - k))
    return r2_base, r2_full, beta[-1], stats.f.sf(f, 1, n - k)


def main():
    human = load_human()
    split = get_split()
    files = sorted((OUT / "surprisal").glob("*.csv"))
    if not files:
        raise SystemExit("No outputs/surprisal/*.csv found: run compute_priming.py first.")
    runs = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    runs = runs.merge(human, on="target").merge(split[["target", "half"]], on="target")
    runs["config"] = runs["model"] + " | " + runs["format"]
    layer_cols = sorted(c for c in runs.columns if c.startswith("dSim_L"))
    disc, conf = runs[runs.half == "discovery"].copy(), runs[runs.half == "confirmation"].copy()
    models = sorted(runs.model.unique())

    # ---------------------------------------------------------------- A. layer profile
    print("=== A. Layer profile (Spearman r with human priming; layer chosen on discovery) ===")
    prof, best_layer = [], {}
    for m in models:
        ref_fmt = "bare" if "bare" in set(runs[runs.model == m].format) else runs[runs.model == m].format.iloc[0]
        for half, d in [("discovery", disc), ("confirmation", conf)]:
            d = d[(d.model == m) & (d.format == ref_fmt)]
            for c in layer_cols:
                if d[c].isna().all():  # this model has fewer layers than the deepest model in the grid
                    continue
                prof.append(dict(model=m, half=half, layer=int(c[-2:]),
                                 r_human=stats.spearmanr(d[c], d["avgZ_prime"])[0],
                                 r_modelDS=stats.spearmanr(d[c], d["dS"])[0]))
    prof = pd.DataFrame(prof)
    prof.to_csv(OUT / "layer_profile.csv", index=False)
    for m in models:
        pm = prof[(prof.model == m)]
        pd_ = pm[pm.half == "discovery"].set_index("layer")
        best_layer[m] = int(pd_["r_human"].abs().idxmax())  # |r|: sign of similarity-priming link is not assumed
        last = pd_.index.max()
        cf = pm[(pm.half == "confirmation") & (pm.layer == best_layer[m])]["r_human"].iloc[0]
        cl = pm[(pm.half == "confirmation") & (pm.layer == last)]["r_human"].iloc[0]
        print(f"{m}: best layer (discovery) = {best_layer[m]} of {last}; confirmation r at best = {cf:.3f}; at last layer = {cl:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=False)
    for m in models:
        for ax, col, title in [(axes[0], "r_human", "Similarity vs human priming"),
                               (axes[1], "r_modelDS", "Similarity vs model's own dS")]:
            pm = prof[(prof.model == m) & (prof.half == "discovery")]
            ax.plot(pm.layer / pm.layer.max(), pm[col], marker="o", ms=3, label=m)
            ax.set_title(title); ax.set_xlabel("relative depth (0 = embeddings, 1 = last layer)")
            ax.set_ylabel("Spearman r (discovery half)"); ax.axhline(0, color="grey", lw=.5)
    axes[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "layer_profile.png", dpi=200); plt.close(fig)
    print("saved outputs/layer_profile.png")

    # similarity covariate = each model's discovery-chosen layer
    for df in (runs, disc, conf):
        df["dSim"] = [row[f"dSim_L{best_layer[row['model']]:02d}"] for _, row in df.iterrows()]
    runs["dSim"] = runs["dSim"].astype(float)
    covars = ["dSim"] + COVARS

    # ---------------------------------------------------------------- B. discovery grid
    rows = []
    for cfg, d in disc.groupby("config"):
        r, lo, hi = spearman_ci(d["dS"], d["avgZ_prime"], n_boot=300)
        rows.append(dict(config=cfg, model=d.model.iloc[0], format=d.format.iloc[0], n=len(d), r_dS=r, ci_lo=lo, ci_hi=hi))
    grid = pd.DataFrame(rows).sort_values("r_dS", ascending=False)
    print(f"\n=== B. Discovery grid ({len(grid)} configs; rule: {SELECTION_RULE}) — top 12 ===")
    print(grid.head(12).round(3).to_string(index=False))
    grid.to_csv(OUT / "discovery_grid.csv", index=False)

    # ---------------------------------------------------------------- C. confirmation
    best = grid.iloc[0]["config"]
    print(f"\nWinner on discovery: {best}")
    summ = []
    for half, df in [("discovery", disc), ("confirmation", conf)]:
        d = df[df.config == best].dropna(subset=["dS", "avgZ_prime"] + covars)
        r, lo, hi = spearman_ci(d["dS"], d["avgZ_prime"])
        r2b, r2f, b, p = incremental(d, covars)
        summ.append(dict(half=half, n=len(d), spearman=r, ci_lo=lo, ci_hi=hi,
                         R2_covariates=r2b, R2_with_dS=r2f, dR2=r2f - r2b, beta_dS=b, p_dS=p))
    s = pd.DataFrame(summ)
    print("\n=== C. Winner: discovery vs confirmation ===")
    print(s.round(4).to_string(index=False))
    print("covariates:", ", ".join(covars), f"(dSim = layer {best_layer[best.split(' | ')[0]]} of the winning model)")
    s.to_csv(OUT / "winner_summary.csv", index=False)

    # ---------------------------------------------------------------- D. frame / shots (descriptive)
    cd = pd.DataFrame([dict(model=d.model.iloc[0], format=d.format.iloc[0],
                            r_confirmation=stats.spearmanr(d["dS"], d["avgZ_prime"])[0])
                       for _, d in conf.groupby("config")])
    print("\n=== D. Confirmation r by model x format (descriptive; winner already fixed) ===")
    wide = cd.pivot(index="format", columns="model", values="r_confirmation")
    print(wide.round(3).to_string())
    wide.to_csv(OUT / "confirmation_by_format.csv")

    # ---------------------------------------------------------------- F. distributions + noise ceiling
    raw_ml = pd.read_csv(ROOT / "data" / "priming_data" / "en_priming.csv")
    se2 = raw_ml["seZ_RT_related"] ** 2 + raw_ml["seZ_RT_unrelated"] ** 2
    rel_h = 1 - se2.mean() / raw_ml["avgZ_prime"].var()
    print("\n=== F. Do the priming distributions match? (winner config, all 1,000 items) ===")
    print(f"Human item-level reliability (from SEs) = {rel_h:.2f}; max attainable r ~ {np.sqrt(max(rel_h, 0)):.2f}")
    w = runs[runs.config == best].dropna(subset=["dS", "avgZ_prime"])
    hz, mz = z(w["avgZ_prime"]), z(w["dS"])
    rows = []
    for name, v in [("human avgZ_prime", w["avgZ_prime"]), ("model dS (nats)", w["dS"])]:
        rows.append(dict(measure=name, mean=v.mean(), sd=v.std(), d=v.mean() / v.std(),
                         skew=stats.skew(v), ex_kurtosis=stats.kurtosis(v), pct_positive=100 * (v > 0).mean()))
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    ks = stats.ks_2samp(hz, mz)
    print(f"Shape after z-scoring: KS D = {ks.statistic:.3f} (p = {ks.pvalue:.3g}); "
          f"item-level Spearman r = {stats.spearmanr(w['dS'], w['avgZ_prime'])[0]:.3f}; "
          f"sign agreement = {100 * ((w['dS'] > 0) == (w['avgZ_prime'] > 0)).mean():.1f}%")
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    bins = np.linspace(-4, 4, 41)
    ax[0].hist(hz, bins, alpha=.55, density=True, label="human"); ax[0].hist(mz, bins, alpha=.55, density=True, label="model")
    ax[0].set_title("Distributions (z-scored)"); ax[0].legend()
    qs = np.linspace(.01, .99, 99)
    ax[1].plot(np.quantile(hz, qs), np.quantile(mz, qs), ".", ms=4); ax[1].plot([-3, 3], [-3, 3], "k--", lw=.7)
    ax[1].set_xlabel("human quantiles"); ax[1].set_ylabel("model quantiles"); ax[1].set_title("QQ plot")
    ax[2].scatter(hz, mz, s=6, alpha=.4); ax[2].set_xlabel("human priming (z)"); ax[2].set_ylabel("model dS (z)")
    ax[2].set_title(f"Item level (Spearman r = {stats.spearmanr(hz, mz)[0]:.2f})")
    fig.suptitle(best, fontsize=9); fig.tight_layout(); fig.savefig(OUT / "distributions.png", dpi=200); plt.close(fig)
    print("saved outputs/distributions.png")

    # ---------------------------------------------------------------- G. external validation on original SPP
    spp_rows = []
    for kind in ["first", "other"]:
        sdir = OUT / f"spp_{kind}"
        items_f = ROOT / "data" / "spp" / f"spp_items_{kind}.csv"
        if not sdir.exists() or not items_f.exists():
            continue
        hum = pd.read_csv(items_f)
        for f in sorted(sdir.glob("*.csv")):
            m = pd.read_csv(f)
            if f"{m.model.iloc[0]} | {m.format.iloc[0]}" not in set(grid.head(8)["config"]):
                continue  # only the 8 configs chosen on ManyLanguages discovery (rule in run_spp_validation.py)
            m = m.merge(hum[["target", "short_priming", "long_priming", "in_ml_set"]], on="target")
            for subset, d in [("all targets", m), ("targets not in ML set", m[~m.in_ml_set])]:
                for soa in ["short", "long"]:
                    spp_rows.append(dict(prime=kind, config=f"{m.model.iloc[0]} | {m.format.iloc[0]}", subset=subset,
                                         soa=soa, n=len(d), r=stats.spearmanr(d["dS"], d[f"{soa}_priming"])[0]))
    print("\n=== G. External validation on original SPP (configs fixed from ML discovery; no re-selection) ===")
    if spp_rows:
        g = pd.DataFrame(spp_rows)
        g.to_csv(OUT / "spp_validation.csv", index=False)
        print(g.pivot_table(index=["config", "prime", "subset"], columns="soa", values="r").round(3).to_string())
    else:
        top = " ".join(sorted(set(grid.head(6)["format"])))
        print("no SPP outputs yet. Run, for the top discovery configs:")
        print(f"  python scripts/compute_priming.py --models <models> --formats {top} \\")
        print("     --items-file data/spp/spp_items_first.csv --out-dir outputs/spp_first   (and _other)")

    # ---------------------------------------------------------------- E. supervised readout ceiling
    piv = runs.pivot_table(index="target", columns="config", values="dS")
    sim = runs.drop_duplicates(["target", "model"]).pivot_table(index="target", columns="model", values="dSim")
    sim.columns = [f"dSim|{c}" for c in sim.columns]
    X = piv.join(sim).join(human.set_index("target")[COVARS]).join(split.set_index("target")["half"])
    X["y"] = human.set_index("target")["avgZ_prime"]
    X = X.dropna()
    feats = [c for c in X.columns if c not in ("half", "y")]
    tr, te = X[X.half == "discovery"], X[X.half == "confirmation"]
    mu, sd = tr[feats].mean(), tr[feats].std()
    ridge = RidgeCV(alphas=np.logspace(-2, 4, 25), cv=KFold(5, shuffle=True, random_state=SEED))
    ridge.fit((tr[feats] - mu) / sd, tr["y"])
    pred = ridge.predict((te[feats] - mu) / sd)
    base = RidgeCV(alphas=np.logspace(-2, 4, 25), cv=KFold(5, shuffle=True, random_state=SEED))
    base_f = [c for c in feats if c in COVARS or c.startswith("dSim|")]
    base.fit((tr[base_f] - mu[base_f]) / sd[base_f], tr["y"])
    pred_b = base.predict((te[base_f] - mu[base_f]) / sd[base_f])
    print("\n=== E. Supervised readout (ceiling; fit on discovery, scored on confirmation) ===")
    print(f"features: {len(feats)} (all dS configs + per-model similarity + covariates); ridge alpha = {ridge.alpha_:.3g}")
    print(f"confirmation r, covariates + similarity only : {stats.pearsonr(pred_b, te['y'])[0]:.3f}")
    print(f"confirmation r, + all dS configs             : {stats.pearsonr(pred, te['y'])[0]:.3f}")


if __name__ == "__main__":
    main()
