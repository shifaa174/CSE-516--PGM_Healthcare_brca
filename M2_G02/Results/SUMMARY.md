# Initial KPI results — M2_G02

Software used for this run: Python 3.12.5, numpy 2.5.3. Runs on Windows, macOS and Linux. Regenerate everything from the `Code` folder with `python run_all.py ../..` (or `bash run_all.sh ../..`).

Model selection rule: the latent dimension d is chosen per learning route by **validation** C-index; test data is never used to choose anything.

## 1. Paper baseline reproduction

Paper methods rerun from the authors' trained SCCA / GN-SCCA models (`03_brca/results/`) and their PCA files, scored with the authors' own survival code (`03_brca/A_survival_prediction/`, lifelines Cox). Largest per-fold difference from the authors' stored results: 0.00 C-index points.

| Method | Setting | Rerun (ours) | Authors' stored results | arXiv Table 4 |
|---|---|---|---|---|
| PCA | genomics | 67.59 ± 7.66 | 67.59 ± 7.66 | 65.76 ± 9.02 |
| PCA | imaging | 58.13 ± 7.97 | 58.13 ± 7.97 | 61.24 ± 6.39 |
| PCA | concat | 66.47 ± 6.86 | 66.47 ± 6.86 | 66.86 ± 10.40 |
| SCCA/HD | genomics | 68.67 ± 9.34 | 68.67 ± 9.34 | 65.95 ± 10.70 |
| SCCA/HD | imaging | 53.15 ± 10.89 | 53.15 ± 10.89 | 54.00 ± 11.23 |
| SCCA/HD | concat | 67.79 ± 8.26 | 67.79 ± 8.26 | 65.71 ± 7.96 |
| SCCA/PD | genomics | 67.79 ± 3.89 | 67.79 ± 3.89 | 68.34 ± 8.40 |
| SCCA/PD | imaging | 56.99 ± 9.89 | 56.99 ± 9.89 | 55.66 ± 11.14 |
| SCCA/PD | concat | 67.58 ± 4.69 | 67.58 ± 4.69 | 67.37 ± 6.58 |
| SCCA/OPD | genomics | 70.03 ± 7.50 | 70.03 ± 7.50 | 67.63 ± 9.43 |
| SCCA/OPD | imaging | 55.86 ± 10.07 | 55.86 ± 10.07 | 54.19 ± 11.23 |
| SCCA/OPD | concat | 68.98 ± 6.16 | 68.98 ± 6.16 | 65.25 ± 5.61 |
| GN-SCCA/HD | genomics | failed | failed | failed |
| GN-SCCA/HD | imaging | failed | failed | failed |
| GN-SCCA/HD | concat | failed | failed | failed |
| GN-SCCA/PD | genomics | 62.60 ± 10.79 | 62.60 ± 10.79 | 64.94 ± 8.98 |
| GN-SCCA/PD | imaging | 58.10 ± 10.46 | 58.10 ± 10.46 | 56.81 ± 13.31 |
| GN-SCCA/PD | concat | 65.47 ± 10.14 | 65.47 ± 10.14 | 66.42 ± 12.08 |
| GN-SCCA/OPD | genomics | 63.94 ± 11.05 | 63.94 ± 11.05 | 65.52 ± 7.67 |
| GN-SCCA/OPD | imaging | 58.86 ± 9.71 | 58.86 ± 9.71 | 57.87 ± 10.28 |
| GN-SCCA/OPD | concat | 66.85 ± 9.47 | 66.85 ± 9.47 | 68.08 ± 9.88 |

GN-SCCA with HD deflation fails in the authors' own results too. The authors' stored results differ from the arXiv preprint's Table 4, which most likely reflects a revision between preprint and journal version; quote the journal table in the report.

## 2. Downstream KPI — our PGM, test C-index (%), mean ± sd over 5 folds

### 2a. On the paper's PCA input: 50 genomic PCs + 50 imaging PCs (like-for-like with the paper's PCA baseline)

| Method | Multimodal p(z\|x,y) | Genomics only p(z\|x) | Histology only p(z\|y) |
|---|---|---|---|
| CCA-MLE, d = 30 | 65.43 ± 6.58 | 67.31 ± 8.61 | 57.34 ± 6.46 |
| EM (our extension), d = 50 | 63.44 ± 8.09 | 69.95 ± 8.68 | 56.93 ± 8.38 |
| Non-PGM baseline + our Cox | 66.08 ± 6.97 | 67.82 ± 7.38 | 57.64 ± 8.49 |

### 2b. On the reduced input: 50 genomic PCs + 215 raw CellProfiler features

| Method | Multimodal p(z\|x,y) | Genomics only p(z\|x) | Histology only p(z\|y) |
|---|---|---|---|
| CCA-MLE, d = 20 | 65.95 ± 10.55 | 67.26 ± 10.50 | 59.23 ± 3.91 |
| EM (our extension), d = 10 | 67.58 ± 7.76 | 72.11 ± 6.17 | 60.30 ± 5.97 |
| Non-PGM baseline + our Cox | 69.78 ± 6.32 | 67.82 ± 7.38 | 62.78 ± 8.45 |

## 3. Learning KPI — log-likelihood per sample and EM convergence

| Input | Route | d | Train logL/n | Test logL/n | Fit time |
|---|---|---|---|---|---|
| paper PCA 50 + 50 | cca-mle | 30 | -132.16 ± 0.38 | -157.10 ± 4.63 | 0.01 s |
| paper PCA 50 + 50 | em | 50 | -131.56 ± 0.37 | -158.68 ± 5.19 | 0.14 s |
| reduced 50 + 215 | cca-mle | 20 | 107.52 ± 0.79 | 93.02 ± 2.62 | 0.02 s |
| reduced 50 + 215 | em | 10 | 338.06 ± 1.88 | 259.20 ± 14.76 | 2.53 s |

- EM on the paper PCA input: converged on all 5 folds (75–104 iterations).
- EM on the reduced input: converged on 2 of 5 folds; on the other 3 it plateaued and oscillated by at most 0.0014% until the 300-iteration cap (positive-definite projection of Psi in the M-step).

## 4. Inference KPI — posterior differential entropy (nats)

| Input | Route | d | Prior | Multimodal | Genomics only | Histology only |
|---|---|---|---|---|---|---|
| paper PCA 50 + 50 | cca-mle | 30 | 42.57 | 28.39 | 33.64 | 33.64 |
| paper PCA 50 + 50 | em | 50 | 70.95 | 51.12 | 59.63 | 58.44 |
| reduced 50 + 215 | cca-mle | 20 | 28.38 | 16.67 | 20.80 | 20.80 |
| reduced 50 + 215 | em | 10 | 14.19 | 0.42 | 2.21 | 7.73 |

Lower than the prior = the data is informative about Z; withholding a modality raises entropy.

## 5. Representation KPI — latent dimension sweeps

### 5a. Paper PCA input

| d | CCA valid logL/n | CCA valid C-index | CCA test C-index | EM valid logL/n | EM valid C-index | EM test C-index |
|---|---|---|---|---|---|---|
| 5 | -149.94 | 59.32 | 51.51 ± 8.44 | -150.85 | 54.55 | 48.85 ± 10.61 |
| 10 | -150.76 | 60.66 | 55.02 ± 7.04 | -151.75 | 61.37 | 54.33 ± 9.34 |
| 20 | -151.94 | 61.70 | 62.58 ± 6.57 | -153.12 | 62.80 | 61.67 ± 6.49 |
| 30 | -152.46 | 62.31 | 65.43 ± 6.58 | -153.70 | 61.53 | 64.81 ± 4.58 |
| 50 | -152.68 | 59.09 | 61.99 ± 4.12 | -153.96 | 63.14 | 63.44 ± 8.09 |

Chosen on validation C-index: CCA-MLE d = 30, EM d = 50. Test C-index is reported only, never used to choose.

### 5b. Reduced input

| d | CCA valid logL/n | CCA valid C-index | CCA test C-index | EM valid logL/n | EM valid C-index | EM test C-index |
|---|---|---|---|---|---|---|
| 2 | 92.39 | 54.09 | 51.73 ± 11.33 | 259.25 | 58.13 | 62.95 ± 6.31 |
| 5 | 92.32 | 65.86 | 65.66 ± 5.85 | 256.82 | 63.74 | 61.85 ± 6.13 |
| 10 | 91.94 | 67.53 | 66.89 ± 4.99 | 252.84 | 65.45 | 67.00 ± 8.23 |
| 20 | 91.12 | 69.12 | 65.95 ± 10.55 | 246.78 | 58.74 | 63.14 ± 5.36 |
| 30 | 90.67 | 64.78 | 64.05 ± 7.58 | 242.96 | 55.85 | 63.28 ± 8.15 |

Chosen on validation C-index: CCA-MLE d = 20, EM d = 10. Test C-index is reported only, never used to choose.

## 6. Findings

1. **The paper's results reproduce.** Rerunning the authors' trained models gives their numbers (best concatenated: SCCA/OPD 68.98 ± 6.16; PCA 66.47 ± 6.86).
2. **On the same input, our PGM does not beat the paper.** CCA-MLE 65.43 ± 6.58 and EM 63.44 ± 8.09 versus PCA 66.47 ± 6.86 and SCCA/OPD 68.98 ± 6.16; all differences are within one standard deviation.
3. **Held-out likelihood cannot choose d.** Under both learning routes, validation log-likelihood keeps falling as d grows, always picking d = 5 on the paper PCA input and always picking d = 2 on the reduced input, which predicts poorly. We therefore choose d by validation C-index.
4. **EM's likelihood gain depends on the input.** At d = 10, EM vs CCA-MLE validation logL/n is 252.84 vs 91.94 on the reduced input, but -151.75 vs -150.76 on the paper's PCA input. The reduced input contains constant and duplicated CellProfiler columns; there the ridge-regularised CCA route is far from the true maximum-likelihood fit and EM corrects it, while on clean PCA features both routes reach the same fit.

## 7. Caveats

- 95 events in 974 patients: fold-to-fold sd of 5–10 points is intrinsic; no difference here is claimed as significant.
- Validation splits hold about 15–25 events each, so the validation-chosen d is itself noisy.
- Across machines, C-index can differ in the second decimal place (numerical libraries); conclusions are unchanged.
- The paper's SCCA / GN-SCCA rows use the authors' trained models; their training code is MATLAB and is not rerun.

## Figures

- `fig1_graphical_model.png` — the implemented PGM and the CI assumption
- `fig2_em_convergence.png` — EM log-likelihood traces, reduced input
- `fig3_cindex_comparison.png` — C-index by method and modality, reduced input
- `fig3b_cindex_paper_input.png` — C-index by method and modality, paper PCA input
- `fig4_latent_dim_sweep.png` — d sweep, reduced input
- `fig4b_latent_dim_sweep_paper_input.png` — d sweep, paper PCA input
- `fig5_posterior_entropy.png` — latent uncertainty by modality availability
- `fig6_paper_comparison.png` — paper baselines (rerun) vs our PGM on the same input
