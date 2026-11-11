#!/usr/bin/env bash
set -e
ROOT="$(cd "${1:?usage: bash run_all.sh <path to repo containing 03_brca/>}" && pwd)"
OUT="$(cd "$(dirname "$0")/.." && pwd)/Results"
cd "$(dirname "$0")/src"

echo "== 1/6  paper baseline: authors' trained models and PCA input =="
python3 paper_baseline.py --root "$ROOT" --out "$OUT"

echo "== 2/6  latent-dimension sweep, paper PCA input =="
python3 sweep_latent_dim.py --root "$ROOT" --variant paperpca --grid 5 10 20 30 50 \
        --em-iters 120 --seed 42 --out "$OUT" --tag paperpca

echo "== 3/6  our PGM on the paper PCA input, d chosen on validation =="
python3 run_experiments.py --root "$ROOT" --variant paperpca --noise block \
        --d-from "$OUT/m2_latent_dim_sweep_paperpca.json" --seed 42 --out "$OUT" --tag paperpca

echo "== 4/6  latent-dimension sweep, reduced input =="
python3 sweep_latent_dim.py --root "$ROOT" --variant pca --grid 2 5 10 20 30 \
        --em-iters 120 --seed 42 --out "$OUT" --tag main

echo "== 5/6  our PGM on the reduced input, d chosen on validation =="
python3 run_experiments.py --root "$ROOT" --variant pca --noise block \
        --d-from "$OUT/m2_latent_dim_sweep_main.json" --seed 42 --out "$OUT" --tag main

echo "== 6/6  figures and summary =="
python3 make_figures.py "$OUT"
python3 write_summary.py "$OUT"

echo "done -- results in $OUT"
