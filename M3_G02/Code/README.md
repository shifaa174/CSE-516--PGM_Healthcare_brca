# M3 Analysis Code

These scripts extend the existing M2 implementation without replacing it. They read the existing M3 JSON outputs, and can run additional robustness experiments using the existing `data.py`, `pcca.py`, and `survival.py`.

## 1. Install dependencies

From the repository root:

```bash
python -m pip install -r M2_G02/Code/requirements.txt pandas tabulate
```

## 2. Aggregate the existing M3 JSON files

```bash
python M3_G02/Code/src/aggregate_results.py \
  --results M3_G02/Results \
  --out M3_G02/Results/analysis
```

Outputs:
- `run_level_metrics.csv`
- `fold_level_metrics.csv`
- `experiment_summary.csv`
- `cindex_summary_percent.csv`
- `paired_tests_exploratory.csv`
- `ANALYSIS_README.md`

The aggregation script does not rerun experiments. It extracts the saved JSON values.

### Statistical interpretation

The five folds are fixed, and seeds 42-46 change EM initialisation only. Therefore:
- Do not describe the five seeds as five independent datasets.
- Standard deviation across folds is not a patient-level confidence interval.
- Paired tests across five folds are exploratory because the sample of folds is small.
- Do not claim statistical significance based only on these exploratory tests.

## 3. Run additional robustness experiments

```bash
python M3_G02/Code/src/run_robustness.py \
  --root . \
  --variant paperpca \
  --folds 0 1 2 3 4 \
  --methods cca-mle em \
  --d-cca 30 --d-em 50 \
  --seed 42 \
  --out M3_G02/Results/robustness
```

This runs three separate analyses:

1. **Missing modality:** randomly hides either genomics or histology for a specified fraction of test patients. Remaining patients retain both modalities.
2. **Observation noise:** adds Gaussian noise to standardised test features at specified levels.
3. **Training size:** trains on different random fractions of the training split while keeping the test split fixed.

Default levels:
- Missing modality: 0%, 10%, 25%, 50%.
- Gaussian noise SD: 0, 0.1, 0.25, 0.5, 1.0.
- Training fraction: 25%, 50%, 75%, 100%.

The missing-modality experiment is a simulated inference-time stress test, not a claim that real-world missingness is random. The training-size experiment is repeated once per fold/seed by default; for stronger uncertainty estimates, rerun with multiple `--seed` values and aggregate separately.

Important: in the training-size experiment, latent dimension is reduced only if necessary to fit the smaller training subset. This is recorded in the CSV.

## 4. Generate plots

```bash
python M3_G02/Code/src/plot_m3_results.py \
  --analysis M3_G02/Results/analysis \
  --robustness M3_G02/Results/robustness \
  --out M3_G02/Results/plots
```

## 5. Suggested execution order

```bash
# A. Analyse existing JSONs first
python M3_G02/Code/src/aggregate_results.py

# B. Run additional robustness tests
python M3_G02/Code/src/run_robustness.py --root .

# C. Generate plots
python M3_G02/Code/src/plot_m3_results.py
```

## 6. Outputs to inspect before writing the report

- `analysis/ANALYSIS_README.md`
- `analysis/cindex_summary_percent.csv`
- `analysis/paired_tests_exploratory.csv`
- `robustness/robustness_summary.csv`
- `plots/*.png`

Check convergence warnings and inspect all output CSVs before interpreting results. These scripts do not manufacture missing results or guarantee that EM converges under every ablation.
