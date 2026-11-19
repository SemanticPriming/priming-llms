#!/bin/bash
# Score the eight English discovery-selected setups (and mGPT) on the 1,000 German pairs.
# Run from the project root: bash scripts/run_german.sh
PY=${PY:-/opt/miniconda3/bin/python}
run() { m=$1; shift; $PY scripts/compute_priming.py --models "$m" --formats "$@" --items-file data/de/de_items.csv --out-dir outputs/de; }
run EleutherAI/pythia-1b ldt
run HuggingFaceTB/SmolLM2-360M ldt newline
run Qwen/Qwen2.5-0.5B-Instruct after_word
run Qwen/Qwen2.5-0.5B ldt pair_begins
run allenai/OLMo-1B-hf after_word and
run ai-forever/mGPT after_word ldt and
