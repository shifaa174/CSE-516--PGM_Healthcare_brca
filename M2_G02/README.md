# M2_G02 — Milestone 2: Technical Foundation

**Project:** Extending a Multimodal Probabilistic Graphical Model for Breast Cancer Survival Prediction
**Course:** CSE516 — Probabilistic Graphical Models · Instructor: Prof. Dhaval Patel
**Track:** Applied AI
**Group G02:** Dhanya Achari (AU2420163) · Parth Jawaharani (AU2420056) · Shifataj Shaikh (AU2420174) · Vishrut Patel (AU2340216)

---

## What exists at M2

The model from M1 is now code, the pipeline runs end to end on all five folds,
both learning routes are implemented, inference produces output, the baseline
runs under identical conditions, and the first KPI numbers are on the board.

```
data  →  pCCA model  →  {CCA-MLE | EM}  →  p(z|·) inference  →  Cox PH  →  C-index
```

## Reproduce every number in the report

```bash
pip3 install -r Code/requirements.txt
cd Code
bash run_all.sh ../..
```

About 6–8 minutes. Writes all JSON, CSV and figures into `Results/`.
Seed 42; see `Code/README.md` for the full seed and configuration list.

## Contents

```
M2_G02/
├── README.md                         this file
├── Report/
│   └── M2_G02_Report.pdf             the M2 report
├── Video/
│   └── M2_G02_Video_Link.txt         walkthrough incl. live pipeline demo
├── Code/
│   ├── README.md                     setup, run commands, seeds, module map
│   ├── requirements.txt
│   ├── run_all.sh                    one-command reproduction
│   └── src/
│       ├── data.py                   dataset pipeline, leakage control
│       ├── pcca.py                   THE GRAPHICAL MODEL + both learning routes + inference
│       ├── survival.py               Cox PH and Harrell's C-index
│       ├── run_experiments.py        experiment driver, all KPIs
│       ├── sweep_latent_dim.py       representation KPI: latent-dimension sweep
│       └── make_figures.py           report figures
├── Data/
│   └── README.md                     provenance, audit, preparation steps
└── Results/
    ├── SUMMARY.md                    KPI tables and findings in readable form
    ├── m2_results_main.json          all per-fold numbers, main configuration
    ├── m2_latent_dim_sweep.json      latent-dimension sweep
    ├── m2_kpi_table.csv              flat KPI table
    └── fig1..fig5 .png               figures
```

## The model, in short

Latent-variable Gaussian PGM (probabilistic CCA):

```
z   ~ N(0, I_d)
x|z ~ N(W_x z + mu_x, Psi_x)      x = genomics   (50 PCs)
y|z ~ N(W_y z + mu_y, Psi_y)      y = histology  (215 PCs)

p(z, x, y) = p(z) p(x|z) p(y|z)           factorisation
X ⊥ Y | Z                                  the conditional independence
```

The absent edge between X and Y is the model's content. It is imposed in code
by keeping the overall noise covariance block-diagonal, so within-modality
correlation lives in `Psi_x` / `Psi_y` and only cross-modality correlation is
routed through Z.

- **Parameter learning:** closed-form MLE via classical CCA (the base paper's
  route, our baseline) and Expectation-Maximization (our extension).
- **Structure learning:** not performed, by design. The latent dimension `d`
  is the one structural hyperparameter.
- **Inference:** exact closed-form Gaussian conditioning for `p(z|x,y)`,
  `p(z|x)`, `p(z|y)` and MAP. A missing modality is marginalised out exactly —
  nothing is imputed.

## Headline initial results

| | CCA-MLE (baseline) | EM (our extension) | PCA + Cox (non-PGM) |
|---|---|---|---|
| Held-out log-lik / sample | 93.75 ± 2.68 | **259.20 ± 14.76** | — |
| C-index, multimodal (%) | 66.92 ± 5.09 | 67.60 ± 7.81 | **69.78 ± 6.34** |
| Fit time | ~0.1 s | ~20 s | — |

EM fits the data far better but gains nothing distinguishable downstream, and
the simple PCA baseline still leads on the multimodal C-index. Full discussion
and the latent-dimension finding are in `Results/SUMMARY.md`.

## Changes since M1

1. The M1 plan was to select `d` on validation by held-out log-likelihood.
   Implemented and tested — the criterion is monotone in `d` and always picks
   the smallest value, which is prognostically poor. We fix `d = 10` and report
   the sweep instead.
2. Experiments run on the reduced (PCA) feature variant rather than the full
   1000/1075 features, because with 585 training patients the full
   within-modality covariance blocks are not identifiable. The full variant is
   M3 work. See `Data/README.md`.
3. Non-PGM PCA baselines were added beyond the M1-declared CCA baseline, so
   the comparison is not only CCA versus EM.

## Reproducibility notes

- Seed 42 globally; EM initialisation uses `seed + fold`.
- Train-only standardisation; validation used only for the `d` sweep; test
  scored once per configuration.
- Configuration and environment are written into every results JSON.
- Verified by a clean-copy run: `bash run_all.sh ../..` reproduces every
  number in `Results/` exactly.
