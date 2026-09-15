# Dataset

**Source:** TCGA-BRCA (The Cancer Genome Atlas, Breast Invasive Carcinoma),
processed data released with the reference work (Subramanian et al., 2024).

**Location in this repository:** `03_brca/data/` (train/validation/test splits in `03_brca/data/splits_journal/`)

## Description
- **Patients:** 974 with both genomic and histology data
- **Genomic (X):** RNA-seq expression, log-transformed; top 1000 genes by coefficient of variation → 1000 features
- **Histology (Y):** CellProfiler features from whole-slide image patches, summarized with 5-bin histograms → 1075 features
- **Outcome:** survival time and event indicator
- **Splits:** 60% train / 15% validation / 25% test, five random splits

## Access
No separate download is needed: clone this repository and use `03_brca/data/`.
Raw TCGA data is available from the GDC Data Portal: https://portal.gdc.cancer.gov/
