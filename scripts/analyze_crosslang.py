"""Data language x prompt language (and model language) summaries. Run from the project root."""
import re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

R = Path(".")
sp = pd.read_csv(R / "outputs/split.csv")
en = pd.read_csv(R / "data/priming_data/en_priming.csv")[["target_word_unique", "avgZ_prime"]].rename(columns={"target_word_unique": "target"})
de = pd.read_csv(R / "data/priming_data/de_priming.csv")
de["en_target"] = de.en_target.astype(str).str.strip()
de = de.merge(sp[["target", "half"]], left_on="en_target", right_on="target", how="left", suffixes=("", "_en"))
de = de[["target_word_unique", "avgZ_prime", "half"]].rename(columns={"target_word_unique": "target"})
en = en.merge(sp[["target", "half"]], on="target", how="left")
HUM = {"en": en, "de": de}


def rho(path, lang):
    d = pd.read_csv(path).merge(HUM[lang], on="target").dropna(subset=["dS", "avgZ_prime"])
    out = {"all": stats.spearmanr(d.dS, d.avgZ_prime)[0]}
    for h in ["discovery", "confirmation"]:
        y = d[d.half == h]
        out[h] = stats.spearmanr(y.dS, y.avgZ_prime)[0]
    return out


def name(m):
    return re.sub(r"[^A-Za-z0-9.-]+", "_", m)


# ---- the eight setups: English vs German prompt, on English vs German data
EIGHT = [("allenai/OLMo-1B-hf", "after_word", "de_after_word"), ("allenai/OLMo-1B-hf", "and", "de_and"),
         ("Qwen/Qwen2.5-0.5B", "ldt", "de_ldt"), ("HuggingFaceTB/SmolLM2-360M", "ldt", "de_ldt"),
         ("Qwen/Qwen2.5-0.5B", "pair_begins", "de_pair_begins"), ("HuggingFaceTB/SmolLM2-360M", "newline", "newline"),
         ("Qwen/Qwen2.5-0.5B-Instruct", "after_word", "de_after_word"), ("EleutherAI/pythia-1b", "ldt", "de_ldt")]
rows = []
for m, f, fd in EIGHT:
    r = dict(model=m.split("/")[-1], frame=f)
    r["EN prompt, EN data"] = rho(R / f"outputs/surprisal/{name(m)}__{f}.csv", "en")
    r["EN prompt, DE data"] = rho(R / f"outputs/de/{name(m)}__{f}.csv", "de")
    base = "outputs/eight_de_prompt" if fd != "newline" else None
    if fd == "newline":  # frame has no words, so the English-prompt runs are the German-prompt runs
        r["DE prompt, EN data"], r["DE prompt, DE data"] = r["EN prompt, EN data"], r["EN prompt, DE data"]
    else:
        r["DE prompt, EN data"] = rho(R / f"{base}/en/{name(m)}__{fd}.csv", "en")
        r["DE prompt, DE data"] = rho(R / f"{base}/de/{name(m)}__{fd}.csv", "de")
    rows.append(r)
cells = ["EN prompt, EN data", "EN prompt, DE data", "DE prompt, EN data", "DE prompt, DE data"]
flat = pd.DataFrame([{**{"model": r["model"], "frame": r["frame"]}, **{f"{c}|{h}": r[c][h] for c in cells for h in ["all", "confirmation"]}} for r in rows])
flat.to_csv(R / "outputs/crosslang_eight.csv", index=False)
print("EIGHT SETUPS (Spearman rho; confirmation half / all 1000)")
for c in cells:
    print(f"{c:22s} conf {flat[c + '|confirmation'].mean():.3f}   all {flat[c + '|all'].mean():.3f}")
print(flat.round(2).to_string(index=False))

# ---- model language: GPT-2 vs German GPT-2 vs mGPT, English vs German frames, English vs German data
rows = []
for m in ["gpt2", "dbmdz/german-gpt2", "ai-forever/mGPT"]:
    for f in ["after_word", "and", "ldt", "de_after_word", "de_and", "de_ldt"]:
        for lang in ["en", "de"]:
            p = R / f"outputs/xlang/{lang}/{name(m)}__{f}.csv"
            if p.exists():
                r = rho(p, lang)
                rows.append(dict(model=m.split("/")[-1], prompt="DE" if f.startswith("de_") else "EN",
                                 frame=f.replace("de_", ""), data=lang.upper(), rho_all=r["all"], rho_conf=r["confirmation"]))
g = pd.DataFrame(rows)
g.to_csv(R / "outputs/crosslang_models.csv", index=False)
print("\nMODEL x PROMPT x DATA (mean over 3 frames; all 1000 items / confirmation half)")
print(g.groupby(["model", "prompt", "data"])[["rho_all", "rho_conf"]].mean().round(3).unstack("data").to_string())
