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
| `src/sweep_latent_dim.py` | Latent-dimension selection |
| `src/run_experiments.py` | Main experiments and ablations |
| `src/run_robustness.py` | Missing modality, noise and training-size experiments |
| `src/aggregate_results.py` | Result tables and exploratory paired tests |
| `src/final_stats.py` | Mean, SD, 95% CI, paired t-test and Wilcoxon |
| `src/probability_quality.py` | NLL, Brier score and calibration |
| `src/make_figures.py`, `src/plot_results.py`, `src/write_summary.py` | Plots, KPI table and summary |
| `run_all.py` | Runs everything above in order |

## Environment

Python 3.10 to 3.13 (tested on macOS and Windows).

```bash
pip install -r M4_G02/Code/requirements.txt
pip install -r M4_G02/Code/requirements-paper.txt
pip install pandas tabulate
```

## Seeds and settings

- Data: 5 fixed folds (585 train / 146 validation / 243 test)
- Seeds 42 to 46 (EM initialisation); CCA-MLE, Cox and the paper baselines are deterministic
- Latent dimension: chosen on validation C-index (CCA-MLE d = 30, EM d = 50)
- Cox penalty 0.1, ridge 1e-2, EM max 300 iterations, tolerance 1e-7

## Reproduce every result

One command, from the repository root (about 6 to 10 minutes):

```bash
python M4_G02/Code/run_all.py .
```

It runs, in order, writing everything to `M4_G02/Results/`:

1. Paper baselines
2. Latent-dimension sweep (paper PCA input)
3. Latent-dimension sweep and PGM (reduced input)
4. Main results, 5 seeds x 5 folds
5. Ablations: diagonal noise, CCA-initialised EM, ridge 1e-3 and 1e-1, Cox penalty 0.01 and 1.0
6. Robustness: missing modality, observation noise, training size
7. Aggregated tables and exploratory paired tests
8. Final statistics: 95% CI, paired t-test, Wilcoxon
9. Probability quality: NLL, Brier score, calibration
10. Figures, KPI table and summary

Any step can be run on its own from `M4_G02/Code/src`, for example:

```bash
python final_stats.py --results ../../Results --out ../../Results/stats
python probability_quality.py --root ../../.. --out ../../Results/probability_quality
```

On another machine, C-index can differ in the second decimal place; conclusions are unchanged.
