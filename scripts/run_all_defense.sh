#!/bin/bash
# Re-run the analysis scripts with COW defence pacts only (ALLIANCE_DEF=defense). Outputs overwrite results/*; to keep the four-treaty-type outputs, copy results/ to results_all_types/ before running.
cd "$(dirname "$0")/.."
export ALLIANCE_DEF=defense
R=Rscript
: > results/rerun_defense_log.txt
for s in 35 36 37 39 40 42 44 45 46 47 48 49 50 51 52 53 54 55 57 58 59 60 61 62 63 64; do
  f=$(ls scripts/${s}_*.R | head -1)
  start=$(date +%s)
  $R "$f" > "/tmp/rerun_$s.out" 2>&1; code=$?
  end=$(date +%s)
  echo "$s $(basename $f) exit=$code secs=$((end-start))" >> results/rerun_defense_log.txt
  if [ $code -ne 0 ]; then tail -5 "/tmp/rerun_$s.out" >> results/rerun_defense_log.txt; fi
done
echo DONE >> results/rerun_defense_log.txt
