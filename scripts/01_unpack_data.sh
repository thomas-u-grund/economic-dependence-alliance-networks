#!/bin/bash
# Unpack the two compressed intermediate files (they exceed GitHub's per-file size limit when uncompressed).
# Run once from the repository root: bash scripts/01_unpack_data.sh
set -e
cd "$(dirname "$0")/.."
[ -f data/tie_trade.csv ] || gunzip -k data/tie_trade.csv.gz
[ -f data/triads_all_years.csv ] || gunzip -k data/triads_all_years.csv.gz
echo "data/tie_trade.csv and data/triads_all_years.csv are in place."
