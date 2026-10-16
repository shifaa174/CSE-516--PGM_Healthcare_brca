# Code — M2_G02

Our latent-variable PGM (two learning routes, exact inference, Cox survival
head) plus a reproduction of the base paper's baselines using the authors' own
code and trained models.

## Environment

Python 3.10+ (tested on macOS / Python 3.12 and Linux / Python 3.13).

**Required** — the whole pipeline runs with only these:

```bash
pip3 install -r requirements.txt          # numpy, scipy, matplotlib
```

**Optional** — lets the paper baseline use the authors' own survival code
(`03_brca/A_survival_prediction/`, lifelines Cox) instead of ours:

```bash
pip3 install -r requirements-paper.txt    # lifelines, scikit-survival, pandas, scikit-learn
```

`paper_baseline.py` detects these automatically. If they are missing it falls
back to our Cox model on the authors' trained models, so the pipeline always
runs. With them installed, it reproduces the authors' stored results exactly.

## Reproduce everything with one command

From the `Code/` folder, on any operating system:

```bash
python run_all.py ../..
```

On macOS / Linux, `python3 run_all.py ../..` or `bash run_all.sh ../..` do the
same. `../..` is the repository root (the folder containing `03_brca/`). Six steps:

1. Paper baseline — the authors' trained SCCA / GN-SCCA models (HD, PD, OPD)
   and their PCA input, all 5 folds
2. Latent-dimension sweep on the paper's PCA input (50 + 50 PCs)
3. Our PGM on the paper's PCA input, d chosen on validation
4. Latent-dimension sweep on the reduced input (50 PCs + 215 raw features)
5. Our PGM on the reduced input, d chosen on validation
6. Figures and `Results/SUMMARY.md`

Roughly 10–20 minutes on a laptop.

## Windows

Works in PowerShell or Command Prompt. Differences from macOS / Linux:

- Use `python` (or `py`) instead of `python3`, and `pip` instead of `pip3`
- Use `python run_all.py ../..`, not `run_all.sh` (bash is not built into Windows)
- Paths with forward slashes, as written above, work as-is

```powershell
cd CSE-516--PGM_Healthcare_brca\M2_G02\Code
pip install -r requirements.txt
pip install -r requirements-paper.txt
python run_all.py ../..
```

If `requirements-paper.txt` fails to install, the pipeline still runs: the
paper baseline falls back to our Cox model.

## Model-selection rule

The latent dimension d is chosen **per learning route by validation C-index**
(`sweep_latent_dim.py` writes `chosen_d`; `run_experiments.py --d-from` reads
it). Test data is never used to choose anything. Held-out likelihood is not
used to choose d, because it falls monotonically as d grows and would always
pick the smallest value.

## Seeds and determinism

| Setting | Value |
|---|---|
| Global seed | `--seed 42` |
| EM initialisation | `seed + fold`, i.e. 42–46 |
| CCA-MLE route | deterministic (closed form) |
| Cox PH | deterministic (L-BFGS from zero) |
| Data splits | fixed, shipped with the base paper |

Same machine, same seed → identical numbers. Across machines, C-index can
differ in the second decimal place; conclusions are unchanged.

## Feature inputs (`--variant`)

| Variant | Genomics X | Histology Y | Used for |
|---|---|---|---|
| `paperpca` | 50 genomic PCs | 50 imaging PCs | like-for-like with the paper's PCA baseline |
| `pca` | 50 genomic PCs | 215 raw CellProfiler features | reduced input |
| `loghist` | 1000 genes | 1075 log-histogram features | the paper's SCCA input (M3) |
| `histogram` | 1000 genes | 1075 histogram features | full features (M3) |

## Module map

| File | Role |
|---|---|
| `src/data.py` | loads the splits for each variant; train-only standardisation |
| `src/pcca.py` | the graphical model: CCA-MLE route, EM route, exact inference, log-likelihood |
| `src/survival.py` | our Cox PH and Harrell's C-index |
| `src/paper_baseline.py` | reruns the paper's methods from the authors' trained models (authors' Cox if available, else ours) |
| `src/sweep_latent_dim.py` | d sweep: validation log-likelihood, validation C-index, test C-index; writes `chosen_d` |
| `src/run_experiments.py` | 5 folds × 2 routes × 3 modality settings + non-PGM baselines |
| `src/make_figures.py` | the report figures |
| `src/write_summary.py` | builds `Results/SUMMARY.md` and `m2_kpi_table.csv` from the JSON results |
| `run_all.py` | runs all six steps; works on Windows, macOS and Linux |
| `run_all.sh` | the same six steps as a bash script (macOS / Linux) |

## Authors' code used

- `03_brca/A_survival_prediction/` — the authors' Python survival pipeline,
  imported unchanged by `paper_baseline.py` when its libraries are installed
- `03_brca/results/*_best.mat` — the authors' trained SCCA / GN-SCCA
  projection matrices (1000 × 50 genomics, 1075 × 50 histology per fold)
- `03_brca/survival_results/*.npz` — the authors' stored C-indices, used only
  as the reproduction check
- `01_methods/` — the authors' MATLAB training code; not rerun (MATLAB)

Source: Subramanian, Syeda-Mahmood and Do, *Artificial Intelligence in
Medicine* 149:102787, 2024; https://github.com/svaishnavi411/modelling_cca_brca_survival

## Fair-comparison guarantees

Within a run, every method sees the same folds, the same train-only
standardisation, the same Cox penalty (0.1) and the same test split. The
validation split is used only to choose d.
