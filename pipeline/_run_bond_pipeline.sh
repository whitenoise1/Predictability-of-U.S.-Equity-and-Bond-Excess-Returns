#!/bin/sh
# Drive the all-maturity bond forecast re-run: 60 (config x sample x maturity)
# slices at 7-way concurrency, then merge+G1, then Hansen SPA. Run from repo root.
set -e
mkdir -p /tmp/bondrun
: > /tmp/bondrun/jobs.txt
for c in C-ENET C-RIDGE C-PLS; do
  for s in full pre_2008 post_2008 post_2016; do
    for m in bond_1y bond_2y bond_3y bond_5y bond_10y; do
      echo "$c $s $m" >> /tmp/bondrun/jobs.txt
    done
  done
done
echo "=== SLICES START $(date +%H:%M:%S) ($(wc -l < /tmp/bondrun/jobs.txt) slices) ==="
cat /tmp/bondrun/jobs.txt | xargs -P 7 -I {} sh -c 'set -- {}; CFG=$1 SMP=$2 MAT=$3 Rscript working_paper/pipeline/rerun_bond_forecasts.R > /tmp/bondrun/slice_$1_$2_$3.log 2>&1; echo "SLICE-DONE $1 $2 $3 rc=$?"'
echo "=== SLICES END $(date +%H:%M:%S) ==="
echo "=== MERGE+G1 ==="
Rscript working_paper/pipeline/merge_bond_forecasts.R > /tmp/bondrun/merge.log 2>&1 && echo "MERGE-OK" || { echo "MERGE-FAIL"; tail -20 /tmp/bondrun/merge.log; exit 1; }
echo "=== BOND SPA ==="
Rscript working_paper/pipeline/run_stage4_bond_spa.R > /tmp/bondrun/spa.log 2>&1 && echo "SPA-OK" || { echo "SPA-FAIL"; tail -20 /tmp/bondrun/spa.log; exit 1; }
echo "=== BOND PIPELINE DONE $(date +%H:%M:%S) ==="
