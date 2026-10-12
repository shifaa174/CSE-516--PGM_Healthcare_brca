# Code — M2_G02

Implementation of the latent-variable PGM, both parameter-learning routes,
exact inference, the survival head and the baselines.

## Environment

Python 3.10+ (developed on 3.11 and 3.12). Three dependencies:

```bash
pip3 install -r requirements.txt
```

`lifelines` is deliberately **not** used — it fails to build in several
environments, and writing the Cox partial likelihood out explicitly in
`src/survival.py` makes the objective auditable.

## Data

The pipeline reads the processed TCGA-BRCA splits that ship inside this
repository under `03_brca/data/splits_journal/`. Nothing needs downloading.
See `../Data/README.md` for provenance and the audit.

## Reproduce everything with one command

From the `Code/` folder:

```bash
bash run_all.sh ../..
```

`../..` is the repository root (the folder containing `03_brca/`), resolved to an
absolute path by the script. This runs
the main experiment, the latent-dimension sweep and the figures, writing
everything into `../Results/`. Takes roughly 6–8 minutes on a laptop.

## Or run the steps individually

```bash
cd src

# main experiment: 5 folds, d=10, seed 42
python3 run_experiments.py --root ../../.. --variant pca --noise block \
        --d 10 --seed 42 --out ../../Results --tag main

# representation KPI: latent-dimension sweep
python3 sweep_latent_dim.py --root ../../.. --grid 2 5 10 20 30 \
        --em-iters 120 --seed 42 --out ../../Results

# figures
python3 make_figures.py ../../Results
```

## Seeds and determinism

| Setting | Value |
|---|---|
| Global seed | `--seed 42` |
| EM initialisation | `seed + fold`, i.e. 42–46 across the five folds |
| CCA-MLE route | deterministic (closed form, no seed needed) |
| Cox PH | deterministic (L-BFGS from a zero start) |
| Data splits | fixed, shipped with the base paper — not resampled |

Re-running with the same seed on the same machine reproduces every number
exactly. Across machines, C-index values can differ in the second decimal
place due to numerical-library differences; conclusions are unchanged. The seed and
the full configuration are written into each results JSON under `config`.

## Command-line options

| Flag | Default | Meaning |
|---|---|---|
| `--root` | required | repository root containing `03_brca/` |
| `--variant` | `pca` | `pca` = 50 genomic + 215 histology PCs; `histogram` = 1000 + 1075 raw features |
| `--noise` | `block` | `block` = full within-modality covariance; `diag` = variances only |
| `--d` | `10` | latent dimension |
| `--select-d` | off | pick `d` on validation by log-likelihood (see note below) |
| `--em-iters` | `300` | EM iteration cap |
| `--em-init` | `random` | `random` or `cca` warm start |
| `--seed` | `42` | global seed |
| `--folds` | `0 1 2 3 4` | which folds to run |

**Note on `--select-d`:** we implemented it, but held-out log-likelihood turns
out to be monotonically decreasing in `d`, so it always returns the smallest
value in the grid, which is prognostically poor. We therefore fix `d = 10` and
report the full sweep instead. This is a finding, not a bug — see
`fig4_latent_dim_sweep.png` and the report.

## Module map

| File | Role | Key contents |
|---|---|---|
| `src/data.py` | dataset pipeline | split loading; `standardise` fits z-scoring on **train only** (the leakage control) |
| `src/pcca.py` | **the graphical model** | `ProbabilisticCCA`: generative model, `fit_cca_mle` (baseline route), `fit_em` (our extension), `posterior` / `map_estimate` / `posterior_entropy` (inference), `score` (evidence) |
| `src/survival.py` | survival head + metric | `CoxPH` (L2-penalised Breslow partial likelihood) and `c_index` (Harrell) |
| `src/run_experiments.py` | experiment driver | all folds × both routes × three modality settings, plus the non-PGM baselines; writes the KPI JSON |
| `src/sweep_latent_dim.py` | representation KPI | sweeps `d`, records validation log-likelihood and test C-index side by side |
| `src/make_figures.py` | figures | the five report figures |

## Where each M2 requirement lives

- **Graphical model** — `pcca.py`, module docstring and `_assemble()`. The
  block-diagonal structure of `Psi` is where `X ⊥ Y | Z` is imposed.
- **Parameter learning** — `fit_cca_mle` (closed-form MLE) and `fit_em`
  (E-step / M-step marked in comments).
- **Structure learning** — not performed, by design. `d` is the one structural
  hyperparameter; see `sweep_latent_dim.py` and the report's justification.
- **Inference** — `posterior()`. A missing modality is marginalised out
  exactly by omitting its term; nothing is imputed.
- **Baseline** — `fit_cca_mle` is the base paper's route; `run_experiments.py`
  additionally runs PCA-concat / PCA-genomics / PCA-histology + Cox under
  identical folds, standardisation, latent dimension and Cox penalty.
- **Initial KPI results** — `../Results/`.

## Fair-comparison guarantees

Every method in a given run sees the same folds, the same train-only
standardisation, the same latent dimension, the same Cox penalty (0.1) and the
same test split. Nothing is tuned on test. The validation split is used only
for the latent-dimension sweep.
