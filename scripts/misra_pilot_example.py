"""BERT surprisal of PILOT after a related vs unrelated cue (Misra, Ettinger, & Rayz, 2020, Table 1 stimuli).

Stimulus: "[CLS] <cue>. I wanted to become a [MASK]. [SEP]"; Surp(T|C) = -log2 P(T | C).
Facilitation F = Surp(T | unrelated) - Surp(T | related).
Run: /opt/miniconda3/bin/python scripts/misra_pilot_example.py
"""
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

for name in ["bert-base-uncased", "bert-large-uncased"]:
    tok = AutoTokenizer.from_pretrained(name)
    model = AutoModelForMaskedLM.from_pretrained(name).eval()
    tid = tok.convert_tokens_to_ids("pilot")
    res = {}
    for cue in ["airplane", "table"]:
        enc = tok(f"{cue}. I wanted to become a {tok.mask_token}.", return_tensors="pt")
        pos = (enc.input_ids[0] == tok.mask_token_id).nonzero().item()
        with torch.no_grad():
            logp = torch.log_softmax(model(**enc).logits[0, pos], -1)
        p = logp[tid].exp().item()
        res[cue] = (p, -logp[tid].item() / torch.log(torch.tensor(2.0)).item())
        rank = (logp > logp[tid]).sum().item() + 1
        print(f"{name:20s} cue={cue:9s} P(pilot)={p:.5f} surprisal={res[cue][1]:.2f} bits rank={rank}")
    print(f"{name:20s} F = {res['table'][1] - res['airplane'][1]:.2f} bits\n")
