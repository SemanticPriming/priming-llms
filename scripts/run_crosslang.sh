#!/bin/bash
# 2x2 check: data language x model/prompt language. English and German versions of three frames,
# scored with English GPT-2, German GPT-2 and the multilingual mGPT on both datasets.
# Run from the project root: bash scripts/run_crosslang.sh
PY=${PY:-/opt/miniconda3/bin/python}
FORMATS="after_word and ldt de_after_word de_and de_ldt"
for m in gpt2 dbmdz/german-gpt2 ai-forever/mGPT; do
  $PY scripts/compute_priming.py --models "$m" --formats $FORMATS --out-dir outputs/xlang/en
  $PY scripts/compute_priming.py --models "$m" --formats $FORMATS --items-file data/de/de_items.csv --out-dir outputs/xlang/de
done
