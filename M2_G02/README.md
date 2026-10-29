# M2_G02 — Milestone 2: Technical Foundation

**Project:** Extending a Multimodal Probabilistic Graphical Model for Breast Cancer Survival Prediction
**Course:** CSE516 — Probabilistic Graphical Models · Instructor: Prof. Dhaval Patel
**Track:** Applied AI
**Group G02:** Dhanya Achari (AU2420163) · Paarth Jawaharani (AU2420056) · Shifataj Shaikh (AU2420174) · Vishrut Patel (AU2340216)

**Base paper:** V. Subramanian, T. Syeda-Mahmood, M. N. Do, "Modelling-based joint embedding of histology and genomics using canonical correlation analysis for breast cancer survival prediction," *Artificial Intelligence in Medicine* 149:102787, 2024. Code and data: https://github.com/svaishnavi411/modelling_cca_brca_survival

---

## What exists at M2

The pipeline runs end to end on all five folds, on the paper's TCGA-BRCA data:

```
paper's data → paper's baselines (authors' trained models, rerun)
             → our pCCA model → {CCA-MLE | EM} → p(z|·) inference → Cox PH → C-index
```

- **Paper baseline:** the paper's PCA, SCCA and GN-SCCA results (HD, PD, OPD
  deflation), rerun from the authors' own trained models and survival code
- **Our model:** the same latent-variable PGM, learned two ways — closed-form
  CCA-MLE and EM (our extension) — with exact posterior and MAP inference
- **Like-for-like comparison:** our model is run on the paper's exact PCA
  input (50 genomic + 50 imaging PCs)

## Reproduce every number

```bash
pip install -r Code/requirements.txt
pip install -r Code/requirements-paper.txt      # optional: authors' survival code
cd Code
python run_all.py ../..
```

Works on Windows, macOS and Linux (on macOS / Linux, `python3` / `pip3` if
`python` is not found). See `Code/README.md` for Windows notes.

Roughly 10–20 minutes. Without the optional libraries the paper baseline falls
back to our Cox model on the authors' trained models, so it always runs.

## Contents

```
M2_G02/
├── README.md
├── Report/        G2_M2_Report_Healthcare_Applied.pdf, G02_M2_Base_Paper_Block_Diagram.pdf
├── Video/         Video.md (link to the video walkthrough)
├── Code/
│   ├── README.md               setup, commands, seeds, module map
│   ├── requirements.txt        numpy, scipy, matplotlib
│   ├── requirements-paper.txt  lifelines, scikit-survival, pandas, scikit-learn (optional)
│   ├── run_all.py              one-command reproduction (any OS)
│   ├── run_all.sh              same, as a bash script
│   └── src/
│       ├── data.py             dataset loading for each input variant
│       ├── pcca.py             the graphical model, both learning routes, inference
│       ├── survival.py         Cox PH and Harrell's C-index
│       ├── paper_baseline.py   paper's methods from the authors' trained models
│       ├── sweep_latent_dim.py latent-dimension sweep and validation choice of d
│       ├── run_experiments.py  5 folds × 2 routes × 3 modality settings + baselines
│       ├── make_figures.py     figures
│       └── write_summary.py    Results/SUMMARY.md and KPI table
├── Data/          README.md: provenance, audit, preparation
└── Results/       SUMMARY.md, JSON, CSV, figures
```

## The model

```
z   ~ N(0, I_d)
x|z ~ N(W_x z + mu_x, Psi_x)      x = genomics
y|z ~ N(W_y z + mu_y, Psi_y)      y = histology

p(z, x, y) = p(z) p(x|z) p(y|z)
X ⊥ Y | Z
```

The conditional independence is imposed by keeping the noise covariance
block-diagonal. Parameter learning: CCA-MLE (baseline route) and EM (our
extension). No structure learning; the latent dimension d is chosen per route
by validation C-index. Inference is exact Gaussian conditioning; a missing
modality is marginalised out, not imputed.

## Results

All numbers, tables and findings: `Results/SUMMARY.md`. Key figure:
`Results/fig6_paper_comparison.png`.

## Changes since M1

1. **Paper baseline rerun from the authors' own code and trained models**, so
   our baseline is the paper's actual method, and checked against the
   authors' stored results.
2. **Our model is also run on the paper's exact PCA input**, so the
   comparison is like-for-like.
3. **d is chosen on validation C-index.** The M1 plan (validation
   log-likelihood) was tested and rejected: likelihood keeps falling as d
   grows and always picks the smallest d.
4. Full 1000 / 1075 features deferred to M3: with 585 training patients the
   full within-modality covariance blocks cannot be estimated.

## Attribution

`01_methods/`, `02_simulations/` and `03_brca/` at the repository root are the
authors' code, trained models and processed data, used unchanged and credited
above. Everything in `M2_G02/Code/` is our implementation.

## Reproducibility

Seed 42; EM initialisation uses seed + fold. Train-only standardisation;
validation used only to choose d; test scored once. Configuration and
environment are written into every results JSON. Across machines, C-index can
differ in the second decimal place; conclusions are unchanged.
