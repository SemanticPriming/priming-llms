"""Compute model-derived priming (surprisal difference) for the English priming items.

For every target, with its related and unrelated cue:

    dS = S(target | unrelated cue) - S(target | related cue),   S = -log P

Positive dS = the model finds the target more predictable after the related
cue, the analogue of the human avgZ_prime (unrelated minus related RT).

Causal LMs (gpt2, pythia, bloom, olmo, ...): teacher-forced log P of the target tokens.
Masked LMs (roberta, bert, ...): PLL-word-l2r (Kauf & Ivanova, 2023) over the target tokens,
with a sentence-final period after the target (MASKED_SUFFIX).

Also saves, for EVERY layer L, the static cosine similarity difference
    dSim_L = cos(cue_related, target) - cos(cue_unrelated, target)
with each word encoded alone (layer 0 = input embeddings), so the analysis can ask
which layer's similarity best matches human priming and the model's own dS.

Format names:
    <frame>        e.g. bare, comma, given_word      (zero-shot)
    <frame>+k<N>   e.g. bare+k5                      (N related example pairs shown first)
Few-shot examples are drawn ONLY from the discovery half, never include the item being
scored, and are identical for the related and unrelated version of an item.

Usage (from the project root):
    python scripts/compute_priming.py --models gpt2 roberta-base --formats bare bare+k5
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoConfig, AutoModel, AutoModelForCausalLM, AutoModelForMaskedLM, AutoTokenizer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, SEED, get_split  # noqa: E402

ITEMS = ROOT / "data" / "priming_data" / "en_priming.csv"
# Masked models are scored with sentence-final punctuation after the target. Without it, a [MASK]
# directly before [SEP] makes BERT put nearly all probability on punctuation, so every target gets a
# huge, nearly constant surprisal (found by checking raw surprisals; see README).
MASKED_SUFFIX = "."
OUT_DIR = ROOT / "outputs" / "surprisal"

# {cue} and {target} are required; the target must come last.
FRAMES = {
    # minimal contexts
    "bare": "{cue} {target}",
    "comma": "{cue}, {target}",
    "newline": "{cue}\n{target}",
    "list": "Words: {cue}, {target}",
    "ldt": "Word 1: {cue}\nWord 2: {target}",
    # relation-implying minimal frames
    "and": "{cue} and {target}",
    "assoc": "{cue} makes me think of {target}",
    # the four pilot carrier sentences, kept to test whether the pilot null was a framing artifact
    "given_word": "Given the word {cue}, the next word is {target}",
    "after_word": "After the word {cue}, comes the word {target}",
    "pair_begins": "The pair begins with {cue}, and continues with {target}",
    "first_word": "If the first word is {cue}, the next word you'll see is {target}",
}


def pick_device():
    if torch.backends.mps.is_available():
        return "mps"
    return "cuda" if torch.cuda.is_available() else "cpu"


def parse_format(name):
    m = re.fullmatch(r"([a-z_]+)(?:\+k(\d+))?", name)
    if not m or m.group(1) not in FRAMES:
        raise SystemExit(f"Unknown format '{name}'. Frames: {list(FRAMES)}; few-shot: <frame>+k<N>")
    return m.group(1), int(m.group(2) or 0)


def build_prefixes(frame, k, items, cue_col, pool):
    """Text that precedes the target for each item, plus the separator before the target."""
    template = FRAMES[frame]
    head = template.split("{target}")[0]
    sep = " " if head.endswith(" ") else ""
    joiner = "\n\n" if "\n" in template else "\n"
    prefixes = []
    for i, row in items.iterrows():
        shots = ""
        if k:
            rng = np.random.default_rng([SEED, i])
            cand = pool[pool["target"] != row["target"]]
            pick = cand.iloc[rng.choice(len(cand), size=k, replace=False)]
            shots = joiner.join(template.format(cue=c, target=t) for c, t in
                                zip(pick["cue_related"], pick["target"])) + joiner
        prefixes.append(shots + head.format(cue=row[cue_col]).rstrip(" "))
    return prefixes, sep


def is_masked(model_name):
    cfg = AutoConfig.from_pretrained(model_name)
    return any(a.endswith("MaskedLM") for a in (cfg.architectures or []))


# ---------------------------------------------------------------- causal
@torch.no_grad()
def causal_surprisal(model, tok, device, prefixes, targets, sep, batch_size=16):
    bos = tok.bos_token_id if tok.bos_token_id is not None else tok.eos_token_id
    seqs, spans = [], []
    for pre, tgt in zip(prefixes, targets):
        ctx = ([bos] if bos is not None else []) + tok(pre, add_special_tokens=False)["input_ids"]
        tid = tok(sep + tgt, add_special_tokens=False)["input_ids"]
        seqs.append(ctx + tid)
        spans.append((len(ctx), len(ctx) + len(tid)))
    out = np.zeros(len(seqs))
    ntok = np.array([b - a for a, b in spans])
    pad = tok.pad_token_id if tok.pad_token_id is not None else 0
    for i in range(0, len(seqs), batch_size):
        chunk = seqs[i:i + batch_size]
        L = max(map(len, chunk))
        ids = torch.full((len(chunk), L), pad)
        att = torch.zeros((len(chunk), L), dtype=torch.long)
        for j, s in enumerate(chunk):
            ids[j, :len(s)] = torch.tensor(s)
            att[j, :len(s)] = 1
        logp = torch.log_softmax(model(input_ids=ids.to(device), attention_mask=att.to(device)).logits.float(), -1).cpu()
        for j, s in enumerate(chunk):
            a, b = spans[i + j]
            out[i + j] = -logp[j, torch.arange(a - 1, b - 1), torch.tensor(s[a:b])].sum().item()
    return out, ntok


# ---------------------------------------------------------------- masked
@torch.no_grad()
def masked_surprisal(model, tok, device, prefixes, targets, sep, batch_size=32):
    """PLL-word-l2r: for target token j, mask tokens j..end of the target word."""
    rows, ntok = [], np.zeros(len(targets), dtype=int)
    suf = tok(MASKED_SUFFIX, add_special_tokens=False)["input_ids"] if MASKED_SUFFIX else []
    for k, (pre, tgt) in enumerate(zip(prefixes, targets)):
        ctx = tok(pre, add_special_tokens=False)["input_ids"]
        tid = tok(sep + tgt, add_special_tokens=False)["input_ids"]
        ntok[k] = len(tid)
        base = [tok.cls_token_id] + ctx + tid + suf + [tok.sep_token_id]
        a = 1 + len(ctx)
        for j in range(len(tid)):
            ids = list(base)
            for m in range(a + j, a + len(tid)):
                ids[m] = tok.mask_token_id
            rows.append((k, ids, a + j, base[a + j]))
    out = np.zeros(len(targets))
    for i in range(0, len(rows), batch_size):
        chunk = rows[i:i + batch_size]
        L = max(len(r[1]) for r in chunk)
        ids = torch.full((len(chunk), L), tok.pad_token_id)
        att = torch.zeros((len(chunk), L), dtype=torch.long)
        for j, r in enumerate(chunk):
            ids[j, :len(r[1])] = torch.tensor(r[1])
            att[j, :len(r[1])] = 1
        logits = model(input_ids=ids.to(device), attention_mask=att.to(device)).logits.float()
        for j, (k, _, pos, true) in enumerate(chunk):
            out[k] += -torch.log_softmax(logits[j, pos], -1)[true].item()
    return out, ntok


# ---------------------------------------------------------------- per-layer similarity
@torch.no_grad()
def layer_embeddings(name, device, words, masked):
    """word -> array [n_layers + 1, hidden]; each word encoded alone, mean over its own tokens.
    A BOS (causal) or CLS/SEP (masked) token is added so the attention-sink position is
    not one of the word's tokens."""
    tok = AutoTokenizer.from_pretrained(name)
    enc = AutoModel.from_pretrained(name).to(device).eval()
    bos = tok.bos_token_id if tok.bos_token_id is not None else tok.eos_token_id
    pad = tok.pad_token_id if tok.pad_token_id is not None else 0
    vecs = {}
    uniq = sorted(set(words))
    for i in range(0, len(uniq), 64):
        chunk = uniq[i:i + 64]
        seqs, spans = [], []
        for w in chunk:
            tid = tok(" " + w, add_special_tokens=False)["input_ids"]
            if masked:
                seqs.append([tok.cls_token_id] + tid + [tok.sep_token_id]); spans.append((1, 1 + len(tid)))
            else:
                pre = [bos] if bos is not None else []
                seqs.append(pre + tid); spans.append((len(pre), len(pre) + len(tid)))
        L = max(map(len, seqs))
        ids = torch.full((len(seqs), L), pad)
        att = torch.zeros((len(seqs), L), dtype=torch.long)
        wmask = torch.zeros((len(seqs), L))
        for j, s in enumerate(seqs):
            ids[j, :len(s)] = torch.tensor(s); att[j, :len(s)] = 1
            wmask[j, spans[j][0]:spans[j][1]] = 1
        hs = enc(input_ids=ids.to(device), attention_mask=att.to(device), output_hidden_states=True).hidden_states
        H = torch.stack(hs, 1).float().cpu()  # [B, layers, L, d]
        m = wmask[:, None, :, None]
        pooled = (H * m).sum(2) / m.sum(2)  # [B, layers, d]
        for w, v in zip(chunk, pooled.numpy()):
            vecs[w] = v
    del enc
    return vecs


def layer_cosines(vecs, a_words, b_words):
    A = np.stack([vecs[w] for w in a_words])  # [N, layers, d]
    B = np.stack([vecs[w] for w in b_words])
    num = (A * B).sum(-1)
    return num / (np.linalg.norm(A, axis=-1) * np.linalg.norm(B, axis=-1) + 1e-12)  # [N, layers]


# ---------------------------------------------------------------- main
def run(model_name, fmt, items, pool, device, vecs, masked, tok, model):
    frame, k = parse_format(fmt)
    fn = masked_surprisal if masked else causal_surprisal
    tg = items["target"].tolist()
    pre_rel, sep = build_prefixes(frame, k, items, "cue_related", pool)
    pre_unrel, _ = build_prefixes(frame, k, items, "cue_unrelated", pool)
    s_rel, ntok = fn(model, tok, device, pre_rel, tg, sep)
    s_unrel, _ = fn(model, tok, device, pre_unrel, tg, sep)

    df = items[["target", "cue_related", "cue_unrelated"]].copy()
    df["model"], df["format"], df["scoring"] = model_name, fmt, "pll_l2r" if masked else "causal"
    df["n_target_tokens"] = ntok
    df["S_related"], df["S_unrelated"] = s_rel, s_unrel
    df["dS"] = s_unrel - s_rel
    dsim = layer_cosines(vecs, items["cue_related"], tg) - layer_cosines(vecs, items["cue_unrelated"], tg)
    for L in range(dsim.shape[1]):
        df[f"dSim_L{L:02d}"] = dsim[:, L]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", required=True)
    ap.add_argument("--formats", nargs="+", default=["bare"])
    ap.add_argument("--limit", type=int, default=None, help="first N items, for smoke tests")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--items-file", default=None,
                    help="alternative item CSV (target, cue_related, cue_unrelated), e.g. data/spp/spp_items_first.csv")
    ap.add_argument("--out-dir", default=None, help="output folder (default outputs/surprisal)")
    args = ap.parse_args()
    for f in args.formats:
        parse_format(f)

    raw = pd.read_csv(ITEMS)
    ml = pd.DataFrame({
        "target": raw["target_word_unique"].astype(str).str.strip(),
        "cue_related": raw["en_cue_word_related"].astype(str).str.strip(),
        "cue_unrelated": raw["en_cue_word_unrelated"].astype(str).str.strip(),
    })
    half = get_split().set_index("target")["half"]
    # few-shot examples always come from the ManyLanguages DISCOVERY half, whatever is being scored
    pool = ml[ml["target"].map(half) == "discovery"].reset_index(drop=True)
    if args.items_file:
        items = pd.read_csv(args.items_file)[["target", "cue_related", "cue_unrelated"]].astype(str)
    else:
        items = ml
    if args.limit:
        items = items.head(args.limit)
    global OUT_DIR
    if args.out_dir:
        OUT_DIR = Path(args.out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    device = pick_device()
    print("device:", device, "| items:", len(items), "| few-shot pool:", len(pool))

    for m in args.models:
        todo = [(f, OUT_DIR / f"{re.sub(r'[^A-Za-z0-9.-]+', '_', m)}__{f}.csv") for f in args.formats]
        todo = [(f, p) for f, p in todo if args.overwrite or not p.exists()]
        if not todo:
            print("skip (all exist):", m)
            continue
        masked = is_masked(m)
        tok = AutoTokenizer.from_pretrained(m)
        if tok.pad_token_id is None and tok.eos_token_id is not None:
            tok.pad_token = tok.eos_token
        model = (AutoModelForMaskedLM if masked else AutoModelForCausalLM).from_pretrained(m).to(device).eval()
        words = pd.concat([items.target, items.cue_related, items.cue_unrelated]).tolist()
        vecs = layer_embeddings(m, device, words, masked)
        for f, path in todo:
            print("running", m, f, flush=True)
            run(m, f, items, pool, device, vecs, masked, tok, model).to_csv(path, index=False)
            print("  wrote", path.name, flush=True)
        del model


if __name__ == "__main__":
    main()
