#!/bin/bash
# Unpack the compressed analysis panel (uncompressed it exceeds GitHub's per-file size limit).
# Run once from the repository root: bash scripts/01_unpack_data.sh
set -e
cd "$(dirname "$0")/.."
[ -f data/tie_trade.csv ] || gunzip -k data/tie_trade.csv.gz
echo "data/tie_trade.csv is in place."
