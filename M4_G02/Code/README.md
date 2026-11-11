# Code — M4_G02

Latent-variable PGM for breast cancer survival prediction: probabilistic CCA
(z -> x, z -> y), learned by CCA-MLE and EM, with exact posterior inference and a
Cox survival model, compared against the base paper's baselines.

## Files

| File | Purpose |
|---|---|
| `src/pcca.py` | Model, learning (CCA-MLE, EM) and inference |
| `src/survival.py` | Cox proportional hazards and C-index |
| `src/data.py` | Data loading and standardisation |
| `src/paper_baseline.py` | Paper baselines (PCA, SCCA, GN-SCCA) |
| `src/run_experiments.py` | Main experiments |
| `src/sweep_latent_dim.py` | Latent-dimension selection |
| `src/run_robustness.py` | Missing modality, noise and training-size experiments |
| `src/aggregate_results.py` | Result tables and paired tests |
| `src/make_figures.py`, `src/plot_m3_results.py`, `src/write_summary.py` | Plots and summaries |
| `run_all.py`, `run_all.sh` | Run the base pipeline end to end |

## Environment

Python 3.10 to 3.13 (tested on macOS and Windows).

```bash
pip install -r M4_G02/Code/requirements.txt
pip install -r M4_G02/Code/requirements-paper.txt
pip install pandas tabulate
```

## Seeds and settings

- Data: 5 fixed folds (585 train / 146 validation / 243 test)
- Seeds 42 to 46 (EM initialisation)
- Latent dimension: CCA-MLE d = 30, EM d = 50 (chosen on validation C-index)
- Cox penalty 0.1, ridge 1e-2, EM max 300 iterations, tolerance 1e-7

## Run

All commands from the repository root.

**1. Base pipeline (baselines, latent-dimension sweep, first results)**
```bash
cd M4_G02/Code && python run_all.py ../.. && cd ../..
```

**2. Main results over 5 seeds** (in bash)
```bash
cd M4_G02/Code/src
bash
python paper_baseline.py --root ../../.. --out ../../Results
for s in 42 43 44 45 46; do
  python run_experiments.py --root ../../.. --variant paperpca --noise block \
    --d-from ../../Results/m2_latent_dim_sweep_paperpca.json \
    --seed $s --out ../../Results --tag paperpca_seed$s
done
exit
cd ../../..
```

**3. Ablations** (in bash)
```bash
cd M4_G02/Code/src
bash
for cfg in "abl_noise_diag:--noise diag" "abl_eminit_cca:--em-init cca" \
           "abl_reg_1e-3:--reg 1e-3" "abl_reg_1e-1:--reg 1e-1" \
           "abl_pen_0.01:--penalizer 0.01" "abl_pen_1.0:--penalizer 1.0"; do
  tag=${cfg%%:*}; flag=${cfg#*:}
  for s in 42 43 44 45 46; do
    python run_experiments.py --root ../../.. --variant paperpca \
      --d-from ../../Results/m2_latent_dim_sweep_paperpca.json \
      --seed $s --out ../../Results --tag ${tag}_seed$s $flag
  done
done
exit
cd ../../..
```

**4. Robustness, tables and plots**
```bash
python M4_G02/Code/src/run_robustness.py --root . --out M4_G02/Results/robustness
python M4_G02/Code/src/aggregate_results.py --results M4_G02/Results --out M4_G02/Results/analysis
python M4_G02/Code/src/plot_m3_results.py --analysis M4_G02/Results/analysis \
  --robustness M4_G02/Results/robustness --out M4_G02/Results/plots
```

**5. Final statistics and probability quality**
```bash
python M4_G02/Code/src/final_stats.py --results M3_G02/Results --out M4_G02/Results/stats
python M4_G02/Code/src/probability_quality.py --root . --out M4_G02/Results/probability_quality
```
`final_stats.py` gives mean, SD and 95% CI per method, and paired t-test and
Wilcoxon tests across the 5 folds (seed results averaged per fold).
`probability_quality.py` gives held-out partial NLL, Brier score at 3 and 5 years,
integrated Brier score and calibration (reliability plot, ECE).

All outputs are written to `M4_G02/Results/`.
