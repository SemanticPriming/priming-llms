#!/bin/bash
# The eight English discovery-selected setups with German-translated frames, on English and German data.
# (SmolLM2 + line break has no words in its frame, so it is unchanged and already scored.)
# Run from the project root: bash scripts/run_eight_de_prompts.sh
PY=${PY:-/opt/miniconda3/bin/python}
run() { m=$1; shift
  $PY scripts/compute_priming.py --models "$m" --formats "$@" --out-dir outputs/eight_de_prompt/en
  $PY scripts/compute_priming.py --models "$m" --formats "$@" --items-file data/de/de_items.csv --out-dir outputs/eight_de_prompt/de; }
run allenai/OLMo-1B-hf de_after_word de_and
run Qwen/Qwen2.5-0.5B de_ldt de_pair_begins
run HuggingFaceTB/SmolLM2-360M de_ldt
run Qwen/Qwen2.5-0.5B-Instruct de_after_word
run EleutherAI/pythia-1b de_ldt
