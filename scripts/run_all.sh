#!/bin/bash
# Re-run the whole analysis in order with the Correlates of War defence-pact definition of alliance.
# Outputs overwrite results/*. Run from anywhere: bash scripts/run_all.sh
cd "$(dirname "$0")/.."
export ALLIANCE_DEF=defense
bash scripts/01_unpack_data.sh   # script 02 rebuilds data/ from raw/ and is not needed after unpacking
: > results/run_all_log.txt
for f in scripts/0[3-9]_*.R scripts/[1-9][0-9]_*.R; do
  start=$(date +%s)
  Rscript "$f" > "/tmp/run_$(basename "$f").out" 2>&1; code=$?
  echo "$(basename "$f") exit=$code secs=$(( $(date +%s) - start ))" >> results/run_all_log.txt
  if [ $code -ne 0 ]; then tail -5 "/tmp/run_$(basename "$f").out" >> results/run_all_log.txt; fi
done
echo DONE >> results/run_all_log.txt
