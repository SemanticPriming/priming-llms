"""Build item files from the original Semantic Priming Project (Hutchison et al., 2013) LDT
item analysis, for use as a SECOND, external validation set.

Each SPP target has a first-associate prime and an 'other' associate prime, each with a
re-paired unrelated prime, and human priming at a short (200 ms) and a long (1200 ms) SOA.
We write one file per prime type with the columns compute_priming.py expects, plus the
human outcomes, and flag targets that also appear in the ManyLanguages English set.

Source file (not copied into the repo): pass its path as the first argument.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
src = Path(sys.argv[1] if len(sys.argv) > 1 else
           ROOT.parents[2] / "2_projects/SPAML/pub_SPAML-PSA/06_Analysis/pre_registered/spp_ldt_prime.xlsx")
d = pd.read_excel(src)
ml = set(pd.read_csv(ROOT / "data/priming_data/en_priming.csv")["target_word_unique"])
out = ROOT / "data" / "spp"
out.mkdir(parents=True, exist_ok=True)
for kind, rel, unrel in [("first", "prime_1st_assoc", "prime_1st_unrelated"),
                         ("other", "prime_other_related", "prime_other_unrelated")]:
    f = pd.DataFrame({
        "target": d["target"].str.lower().str.strip(),
        "cue_related": d[rel].str.lower().str.strip(),
        "cue_unrelated": d[unrel].str.lower().str.strip(),
        "short_priming": d[f"short_{'first' if kind == 'first' else 'other'}_priming"],
        "long_priming": d[f"long_{'first' if kind == 'first' else 'other'}_priming"],
    })
    f["in_ml_set"] = f["target"].isin(ml)
    f.to_csv(out / f"spp_items_{kind}.csv", index=False)
    print(kind, len(f), "items;", int(f.in_ml_set.sum()), "targets also in ML set")
