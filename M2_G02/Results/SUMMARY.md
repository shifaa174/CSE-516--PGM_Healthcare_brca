# Initial KPI results — M2_G02

Configuration: reduced (PCA) feature variant, 50 genomic + 215 histology PCs,
latent dimension d = 10, full block-diagonal noise, Cox penalty 0.1, seed 42,
5 folds. Regenerate with `cd ../Code && bash run_all.sh ../..`.

Raw numbers: `m2_results_main.json`, `m2_latent_dim_sweep.json`,
`m2_kpi_table.csv`.

## Downstream KPI — test C-index (%), mean ± sd over 5 folds

| Method | Multimodal p(z\|x,y) | Genomics only p(z\|x) | Histology only p(z\|y) |
|---|---|---|---|
| CCA-MLE (baseline route) | 66.92 ± 5.09 | 72.11 ± 6.75 | 57.88 ± 6.92 |
| EM (our extension) | 67.60 ± 7.81 | 72.05 ± 6.22 | 60.27 ± 5.98 |
| PCA + Cox (non-PGM baseline) | 69.78 ± 6.34 | 67.81 ± 7.33 | 62.72 ± 8.51 |

## Learning KPI — log-likelihood per sample

| Route | Train | Held-out test | Iterations | Fit time |
|---|---|---|---|---|
| CCA-MLE | 106.3 ± 0.8 | 93.75 ± 2.68 | 1 (closed form) | ~0.1 s |
| EM | 338.1 ± 1.9 | 259.20 ± 14.76 | 103–300 | ~20 s |

EM's log-likelihood increased monotonically on every fold.

## Inference KPI — posterior differential entropy (nats)

Prior entropy at d = 10 is 14.19. Lower = more informative posterior.

| Route | Multimodal | Genomics only | Histology only |
|---|---|---|---|
| CCA-MLE | 6.99 | 9.36 | 9.36 |
| EM | 0.55 | 2.06 | 8.01 |

Withholding either modality raises latent uncertainty, as expected.

## Representation KPI — latent dimension sweep

| d | CCA-MLE valid logL/n | CCA-MLE C-index | EM valid logL/n | EM C-index |
|---|---|---|---|---|
| 2 | 92.39 | 51.71 ± 11.34 | 259.25 | 62.96 ± 6.34 |
| 5 | 92.32 | 65.66 ± 5.76 | 256.82 | 61.78 ± 6.22 |
| 10 | 91.94 | 66.92 ± 5.09 | 252.84 | 67.02 ± 8.31 |
| 20 | 91.12 | 65.91 ± 10.45 | 246.78 | 63.20 ± 5.30 |
| 30 | 90.67 | 64.09 ± 7.50 | 242.96 | 63.31 ± 8.16 |

## Three findings to carry into the report

**1. EM fits the data far better than the closed-form route.** Held-out
log-likelihood 259.2 versus 93.8 nats per sample, and monotone convergence on
every fold. On the learning KPI the extension clearly wins.

**2. The better fit does not buy much prognostic value.** Multimodal C-index
is 67.60 versus 66.92 — well inside one standard deviation, so no difference
can be claimed. A large likelihood gain produced essentially no change in
downstream ranking. This is the M1 hypothesis about likelihood and prognosis
diverging, now with numbers behind it.

**3. Held-out likelihood does not identify a useful latent dimension.**
Validation log-likelihood decreases monotonically in d under both routes, so
likelihood-based selection always returns the smallest d in the grid — which
gives the *worst* C-index for the CCA route (51.7 at d = 2). C-index instead
peaks near d = 10. The unsupervised model-selection criterion we planned at M1
is therefore unusable on its own; we fix d = 10 and report the sweep. This is
a methodological finding, and it is the single most interesting result at M2.

## Honest caveats

- Genomics alone outperforms the multimodal posterior under both routes. The
  second modality is not currently adding prognostic signal.
- The non-PGM PCA baseline beats both PGM routes in the multimodal setting
  (69.78). We do not yet beat the simple baseline, and say so.
- With 95 events in 974 patients, fold-to-fold standard deviations of 5–10
  points are intrinsic. Most differences here are not statistically
  distinguishable, and none is claimed as significant.
- Results are on the reduced feature variant; the full 1000/1075 variant is
  M3 work (see `../Data/README.md` for why).

## Figures

| File | Shows |
|---|---|
| `fig1_graphical_model.png` | the implemented PGM, observed vs latent, the CI assumption |
| `fig2_em_convergence.png` | EM log-likelihood traces, all folds |
| `fig3_cindex_comparison.png` | C-index by method × modality setting, with baselines |
| `fig4_latent_dim_sweep.png` | likelihood vs C-index as d varies — finding 3 |
| `fig5_posterior_entropy.png` | latent uncertainty by modality availability |
