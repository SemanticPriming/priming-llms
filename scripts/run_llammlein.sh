#!/bin/bash
# German-only LLaMmlein with the German frames, on German and English words.
# A watchdog stops each run once its outputs exist, in case the process hangs on exit.
# Usage: M=LSX-UniWue/LLaMmlein_120M bash scripts/run_llammlein.sh   (default is the 120M model)
#
# NOTE FOR THE PAPER: the talk used the 120M model because the 1B download (4.4 GB) was too slow
# on the conference-prep connection. The English comparison model is OLMo-1B, so the sizes do not
# match. For the paper, run LSX-UniWue/LLaMmlein_1B (same size as OLMo-1B and Pythia-1B), and
# consider other German or multilingual models (bloom-1b5-clp-german, gbert-base, EuroLLM-1.7B).
PY=${PY:-/opt/miniconda3/bin/python}
M=${M:-LSX-UniWue/LLaMmlein_120M}
TAG=$(basename "$M")
FR="de_after_word de_and de_ldt"
run() { out=$1; shift
  $PY scripts/compute_priming.py --models $M --formats $FR --out-dir "$out" "$@" &
  pid=$!
  while kill -0 $pid 2>/dev/null; do
    [ "$(ls "$out" 2>/dev/null | grep -c LLaMmlein)" -ge 3 ] && { sleep 5; kill $pid 2>/dev/null; break; }
    sleep 10
  done; }
run outputs/llammlein/$TAG/de --items-file data/de/de_items.csv
run outputs/llammlein/$TAG/en
