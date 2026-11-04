"""LoRA fine-tune a causal LM so its priming contrast dS predicts human item-level priming.

SUPERVISED CEILING, not the model's own priming: the adapter is trained on human avgZ_prime, so
a good score shows how alignable the model is, not what it does by itself (cf. the ridge readout
in analyze.py).

    dS = log P(target | related cue) - log P(target | unrelated cue)      (nats, differentiable)
    prediction = a * dS + b                                               (a, b learned, init by OLS)
    loss = MSE(prediction, z-scored avgZ_prime)

Data use (nothing is chosen on confirmation or SPP):
    train  = 80% of the ManyLanguages DISCOVERY half
    val    = other 20% of discovery (picks the best epoch)
    test   = ManyLanguages CONFIRMATION half (primary)
    extern = original SPP, first and other associates, short and long SOA

Control: --shuffle-labels permutes the human labels within the training items. A real gain must beat
this (and the untrained model), otherwise it is capacity or noise.

Usage (project root):
    python scripts/lora_align.py --model EleutherAI/pythia-410m
    python scripts/lora_align.py --model EleutherAI/pythia-410m --shuffle-labels
"""
import argparse
import copy
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from peft import LoraConfig, get_peft_model
from scipy import stats
from transformers import AutoModelForCausalLM, AutoTokenizer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, SEED, get_split, load_human  # noqa: E402
from compute_priming import FRAMES, ITEMS, build_prefixes, parse_format, pick_device  # noqa: E402

OUT_DIR = ROOT / "outputs" / "lora"


def encode(tok, prefixes, targets, sep):
    bos = tok.bos_token_id if tok.bos_token_id is not None else tok.eos_token_id
    seqs, spans = [], []
    for pre, tgt in zip(prefixes, targets):
        ctx = ([bos] if bos is not None else []) + tok(pre, add_special_tokens=False)["input_ids"]
        tid = tok(sep + tgt, add_special_tokens=False)["input_ids"]
        seqs.append(ctx + tid)
        spans.append((len(ctx), len(ctx) + len(tid)))
    return seqs, spans


def target_logprob(model, tok, device, seqs, spans):
    """Sum of log P over the target tokens, one value per sequence (keeps the graph)."""
    pad = tok.pad_token_id if tok.pad_token_id is not None else 0
    L = max(map(len, seqs))
    ids = torch.full((len(seqs), L), pad)
    att = torch.zeros((len(seqs), L), dtype=torch.long)
    for j, s in enumerate(seqs):
        ids[j, :len(s)] = torch.tensor(s)
        att[j, :len(s)] = 1
    logp = torch.log_softmax(model(input_ids=ids.to(device), attention_mask=att.to(device)).logits.float(), -1)
    out = []
    for j, (a, b) in enumerate(spans):
        tgt = torch.tensor(seqs[j][a:b], device=device)
        out.append(logp[j, torch.arange(a - 1, b - 1, device=device), tgt].sum())
    return torch.stack(out)


class Scorer:
    """Pre-tokenises related/unrelated versions of a set of items and returns dS in batches."""

    def __init__(self, tok, items, frame, pool, device):
        tg = items["target"].tolist()
        pre_r, sep = build_prefixes(frame, 0, items, "cue_related", pool)
        pre_u, _ = build_prefixes(frame, 0, items, "cue_unrelated", pool)
        self.rel = encode(tok, pre_r, tg, sep)
        self.unrel = encode(tok, pre_u, tg, sep)
        self.tok, self.device, self.n = tok, device, len(items)

    def dS(self, model, idx):
        r = target_logprob(model, self.tok, self.device, [self.rel[0][i] for i in idx], [self.rel[1][i] for i in idx])
        u = target_logprob(model, self.tok, self.device, [self.unrel[0][i] for i in idx], [self.unrel[1][i] for i in idx])
        return r - u

    @torch.no_grad()
    def all_dS(self, model, bs=32):
        model.eval()
        out = []
        for i in range(0, self.n, bs):
            out.append(self.dS(model, list(range(i, min(i + bs, self.n)))).cpu())
        return torch.cat(out).numpy()


def rho(x, y):
    return stats.spearmanr(x, y)[0]


def best_format_for(model_name):
    g = pd.read_csv(ROOT / "outputs" / "discovery_grid.csv")
    g = g[(g.model == model_name) & ~g.format.str.contains(r"\+k")]
    return g.iloc[0]["format"], g.iloc[0]["r_dS"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="EleutherAI/pythia-410m")
    ap.add_argument("--format", default=None, help="default: model's best discovery format (zero-shot frames only)")
    ap.add_argument("--rank", type=int, default=8)
    ap.add_argument("--epochs", type=int, default=6)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--val-frac", type=float, default=0.2, help="share of the discovery half held out to pick the epoch")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--shuffle-labels", action="store_true")
    ap.add_argument("--limit", type=int, default=None, help="use only the first N items per set (smoke test)")
    args = ap.parse_args()

    fmt, disc_r = (args.format, None) if args.format else best_format_for(args.model)
    parse_format(fmt)
    tag = f"{re.sub(r'[^A-Za-z0-9.-]+', '_', args.model)}__{fmt}__r{args.rank}__e{args.epochs}__lr{args.lr:g}__v{args.val_frac:g}__s{args.seed}" + ("__SHUFFLED" if args.shuffle_labels else "")
    print(f"model={args.model} format={fmt} (base discovery r={disc_r}) shuffle={args.shuffle_labels} -> {tag}", flush=True)
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    device = pick_device()

    raw = pd.read_csv(ITEMS)
    ml = pd.DataFrame({"target": raw["target_word_unique"].astype(str).str.strip(),
                       "cue_related": raw["en_cue_word_related"].astype(str).str.strip(),
                       "cue_unrelated": raw["en_cue_word_unrelated"].astype(str).str.strip()})
    ml = ml.merge(load_human()[["target", "avgZ_prime"]], on="target").merge(get_split()[["target", "half"]], on="target")
    disc = ml[ml.half == "discovery"].sample(frac=1, random_state=args.seed).reset_index(drop=True)
    n_val = int(len(disc) * args.val_frac)
    val, train = disc.iloc[:n_val].reset_index(drop=True), disc.iloc[n_val:].reset_index(drop=True)
    test = ml[ml.half == "confirmation"].reset_index(drop=True)
    pool = disc[["target", "cue_related"]]  # unused (zero-shot) but required by build_prefixes
    sets = {"train": train, "val": val, "confirmation": test}
    for kind in ["first", "other"]:
        sets[f"spp_{kind}"] = pd.read_csv(ROOT / "data" / "spp" / f"spp_items_{kind}.csv")
    if args.limit:
        sets = {k: v.head(args.limit).reset_index(drop=True) for k, v in sets.items()}
        train, val, test = sets["train"], sets["val"], sets["confirmation"]

    tok = AutoTokenizer.from_pretrained(args.model)
    if tok.pad_token_id is None:
        tok.pad_token = tok.eos_token
    base = AutoModelForCausalLM.from_pretrained(args.model).to(device)
    scorers = {k: Scorer(tok, v, fmt, pool, device) for k, v in sets.items()}

    # ---- untrained baseline
    base_dS = {k: scorers[k].all_dS(base) for k in scorers}
    ytr = train["avgZ_prime"].to_numpy().copy()
    mu, sd = ytr.mean(), ytr.std()
    if args.shuffle_labels:
        ytr = rng.permutation(ytr)
    ytr_z = (ytr - mu) / sd
    print("baseline Spearman  val=%.3f  confirmation=%.3f" % (rho(base_dS["val"], val.avgZ_prime), rho(base_dS["confirmation"], test.avgZ_prime)), flush=True)

    # ---- affine head initialised by OLS on the untrained model
    x0 = base_dS["train"]
    slope, intercept = np.polyfit(x0, ytr_z, 1)
    a = torch.nn.Parameter(torch.tensor(float(slope), device=device))
    b = torch.nn.Parameter(torch.tensor(float(intercept), device=device))

    cfg = LoraConfig(r=args.rank, lora_alpha=2 * args.rank, lora_dropout=0.05, target_modules="all-linear", task_type="CAUSAL_LM")
    model = get_peft_model(base, cfg)
    model.print_trainable_parameters()
    opt = torch.optim.AdamW([{"params": [p for p in model.parameters() if p.requires_grad], "lr": args.lr, "weight_decay": 0.01},
                             {"params": [a, b], "lr": 1e-2, "weight_decay": 0.0}])
    ytr_t = torch.tensor(ytr_z, dtype=torch.float32, device=device)

    best, best_state, hist = -np.inf, None, []
    for ep in range(1, args.epochs + 1):
        model.train()
        order = rng.permutation(len(train))
        losses = []
        for i in range(0, len(order), args.batch):
            idx = order[i:i + args.batch].tolist()
            pred = a * scorers["train"].dS(model, idx) + b
            loss = torch.nn.functional.mse_loss(pred, ytr_t[idx])
            opt.zero_grad()
            loss.backward()
            opt.step()
            losses.append(loss.item())
        v = rho(scorers["val"].all_dS(model), val.avgZ_prime)
        hist.append((ep, float(np.mean(losses)), v))
        print(f"epoch {ep}: train loss {np.mean(losses):.3f}  val Spearman {v:.3f}", flush=True)
        if v > best:
            best, best_ep = v, ep
            best_state = {k: t.detach().cpu().clone() for k, t in model.state_dict().items() if "lora" in k}
    model.load_state_dict({**model.state_dict(), **best_state})
    print(f"best epoch {best_ep} (val Spearman {best:.3f})", flush=True)

    # ---- final scores
    tuned = {k: scorers[k].all_dS(model) for k in scorers}
    rows = []
    for k, df in sets.items():
        ycol = {"train": "avgZ_prime", "val": "avgZ_prime", "confirmation": "avgZ_prime"}.get(k)
        outcomes = ({"human": df["avgZ_prime"]} if ycol else {"short_SOA": df["short_priming"], "long_SOA": df["long_priming"]})
        for name, y in outcomes.items():
            rows.append(dict(set=k, outcome=name, n=len(df), r_untrained=rho(base_dS[k], y), r_lora=rho(tuned[k], y)))
    res = pd.DataFrame(rows)
    res["gain"] = res["r_lora"] - res["r_untrained"]
    print("\n" + res.round(3).to_string(index=False))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res.assign(model=args.model, format=fmt, rank=args.rank, epochs=args.epochs, lr=args.lr, val_frac=args.val_frac, seed=args.seed, shuffled=args.shuffle_labels,
               best_epoch=best_ep).to_csv(OUT_DIR / f"{tag}__summary.csv", index=False)
    pd.concat([df[["target"]].assign(set=k, dS_untrained=base_dS[k], dS_lora=tuned[k]) for k, df in sets.items()]
              ).to_csv(OUT_DIR / f"{tag}__items.csv", index=False)
    print("saved", tag)


if __name__ == "__main__":
    main()
