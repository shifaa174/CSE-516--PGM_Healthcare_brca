# Data — M2_G02

## Dataset

**Processed TCGA-BRCA multimodal cohort**, as released with the base paper's
implementation repository.

- Base paper: V. Subramanian, T. Syeda-Mahmood, M. N. Do, *Artificial
  Intelligence in Medicine* 149:102787, 2024. doi:10.1016/j.artmed.2024.102787
- Source repository: https://github.com/svaishnavi411/modelling_cca_brca_survival
- Primary source: TCGA Breast Invasive Carcinoma, NCI Genomic Data Commons —
  https://portal.gdc.cancer.gov

**Access:** public, no data-use agreement. The processed splits are already
present in this repository under `03_brca/data/splits_journal/`, so the
pipeline needs no download step.

## Layout used by the pipeline

```
03_brca/data/splits_journal/
  {train,valid,test}_histogram/1000/{0,1,2,3,4}/    full features
  {train,valid,test}_pca/100/{0,1,2,3,4}/           reduced features
      genes.csv     genomic feature matrix        X
      image.csv     histology feature matrix      Y
      times.csv     survival time in days         T
      events.csv    event indicator (1 = death)   delta
      patients.csv  TCGA patient barcodes
      gene_names.csv  gene symbol | Entrez ID
```

## Audit — verified by reading the files

| Property | Value |
|---|---|
| Patients total | **974** |
| Train / validation / test | 585 / 146 / 243 (identical in all 5 folds) |
| Genomic features | **1000** (top-variance genes from RNA-seq, named) |
| Histology features | **1075** (CellProfiler descriptors, 5-bin histograms) |
| Reduced variant | 50 genomic PCs, 215 histology PCs |
| Events (deaths) | **95 of 974 = 9.8%** |
| Censoring | 90.2% right-censored |
| Missing values | none — the released matrices are complete |

Per-split event counts, fold 0: train 52/585 (8.9%), validation 19/146
(13.0%), test 24/243 (9.9%).

Reproduce this audit:

```bash
cd Code/src && python3 -c "
from data import load_fold
d = load_fold('../../..', 0, 'histogram')
for k, v in d.items():
    print(k, v['X'].shape, v['Y'].shape, 'events', int(v['events'].sum()))
"
```

**Why the event rate matters.** At 95 events, the concordance index is the
only sensible prediction metric (accuracy is undefined under 90% censoring),
5-fold reporting with a standard deviation is mandatory, and differences
smaller than one standard deviation cannot be claimed.

## Feature inputs used

| Input | Genomics | Histology | Source files | Used for |
|---|---|---|---|---|
| Paper PCA (`paperpca`) | 50 genomic PCs | 50 imaging PCs | `*_pca/100/k/genes.csv`, `*_log_histogram/1000/k/PCA_image.csv` | like-for-like with the paper's PCA baseline |
| Reduced (`pca`) | 50 genomic PCs | 215 raw CellProfiler features | `*_pca/100/k/genes.csv`, `image.csv` | our earlier M2 run |
| Log-histogram (`loghist`) | 1000 genes | 1075 log-histogram features | `*_log_histogram/1000/k/` | input of the paper's SCCA models |

The 215 raw CellProfiler features include constant and duplicated columns;
this is why EM's likelihood advantage over CCA-MLE appears on the reduced
input but not on the paper's PCA input (see `Results/SUMMARY.md`).

Full 1000 / 1075 features are not used to fit our model at M2: with 585
training patients, the full within-modality covariance blocks are not
identifiable. Planned for M3.

## Authors' trained models

`03_brca/results/{sparse,gnscca}_{hd,pd,opd}_1000_{fold}_log_histogram_best.mat`
hold the authors' trained SCCA and GN-SCCA projection matrices (`u`: 1000 × 50,
`v`: 1075 × 50). `Code/src/paper_baseline.py` applies them to the log-histogram
input to rebuild the paper's embeddings. `03_brca/survival_results/*.npz` hold
the authors' stored C-indices, used as the reproduction check.

## Preparation steps performed

1. Load `genes.csv` and `image.csv` per split and fold.
2. Z-score both modalities using **training-split statistics only**
   (`data.standardise`). Validation and test are transformed with the training
   mean and standard deviation — never their own.
3. No imputation: the released matrices have no missing values. Missing
   *modalities* are handled at inference time by marginalisation, not
   imputation.
4. Patients are never shared across splits; the splits ship that way.

## Licensing and ethics

TCGA data is de-identified and publicly released for research use. Only
already-processed feature matrices are used; no patient-identifiable
information is present. Attribution follows NCI GDC guidelines.
