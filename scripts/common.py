"""Shared data loading and the discovery/confirmation split.

The split lives in outputs/split.csv and is created once (seeded, stratified on
word-frequency tercile). compute_priming.py draws few-shot examples ONLY from the
discovery half, so it needs the split before any model is run.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED = 42
SPLIT_FILE = ROOT / "outputs" / "split.csv"


def load_human():
    h = pd.read_csv(ROOT / "data" / "priming_data" / "en_priming.csv")
    h = h.rename(columns={"target_word_unique": "target"})[["target", "avgZ_prime"]]
    freq = pd.read_csv(ROOT / "data" / "subs_count" / "en" / "subs_DF_phono.csv")
    freq = freq.rename(columns={"unigram": "target"})[["target", "unigram_freq"]]
    freq = freq.drop_duplicates("target")  # phonemizer left repeated unigrams (38 rows)
    semd = pd.read_csv(ROOT / "data" / "semdiv_data" / "en" / "en_semantic_diversity.csv")
    semd = semd.rename(columns={"term": "target"})[["target", "semd"]]
    h = h.merge(freq, on="target", how="left").merge(semd, on="target", how="left")
    h["log_freq"] = np.log10(h["unigram_freq"].clip(lower=1))
    h["length"] = h["target"].str.len()
    return h


def get_split():
    """Return target -> half. Reuses the saved split so every run is identical."""
    if SPLIT_FILE.exists():
        return pd.read_csv(SPLIT_FILE)
    h = load_human().copy()
    rng = np.random.default_rng(SEED)
    h["freq_band"] = pd.qcut(h["log_freq"], 3, labels=False, duplicates="drop")
    h["half"] = ""
    for _, idx in h.groupby("freq_band").groups.items():
        idx = rng.permutation(np.array(idx))
        k = len(idx) // 2
        h.loc[idx[:k], "half"] = "discovery"
        h.loc[idx[k:], "half"] = "confirmation"
    out = h[["target", "freq_band", "half"]]
    SPLIT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(SPLIT_FILE, index=False)
    return out
