#!/bin/bash
# mGPT on the German data. mGPT hangs on exit on this machine after writing its files,
# so a watchdog stops it once all six outputs exist.
PY=${PY:-/opt/miniconda3/bin/python}
OUT=outputs/xlang/de
$PY scripts/compute_priming.py --models ai-forever/mGPT --formats after_word and ldt de_after_word de_and de_ldt \
  --items-file data/de/de_items.csv --out-dir $OUT &
pid=$!
while kill -0 $pid 2>/dev/null; do
  [ "$(ls $OUT | grep -c mGPT)" -ge 6 ] && { sleep 5; kill $pid 2>/dev/null; break; }
  sleep 10
done
