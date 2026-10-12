#!/usr/bin/env bash
# Regenerates every number and figure in the M2 report.
#   usage:  bash run_all.sh /path/to/repo-containing-03_brca
set -e
ROOT="$(cd "${1:?usage: bash run_all.sh <path to repo containing 03_brca/>}" && pwd)"  # absolute
OUT="$(cd "$(dirname "$0")/.." && pwd)/Results"
cd "$(dirname "$0")/src"

echo "== 1/3  main experiment (5 folds, d=10, seed=42) =="
python3 run_experiments.py --root "$ROOT" --variant pca --noise block \
        --d 10 --seed 42 --out "$OUT" --tag main

echo "== 2/3  latent-dimension sweep =="
python3 sweep_latent_dim.py --root "$ROOT" --grid 2 5 10 20 30 \
        --em-iters 120 --seed 42 --out "$OUT"

echo "== 3/3  figures =="
python3 make_figures.py "$OUT"

echo "done -- results in $OUT"
