# M3_G02: Milestone 3, Results and Analysis

**Project:** Extending a Multimodal Probabilistic Graphical Model for Breast Cancer Survival Prediction
**Course:** CSE516 Probabilistic Graphical Models, Ahmedabad University · Instructor: Prof. Dhaval Patel
**Track:** Applied AI (Healthcare)

## Team
| Name | Enrollment No. |
|---|---|
| Dhanya Achari | AU2420163 |
| Shifataj Shaikh | AU2420174 |
| Paarth Jawaharani | AU2420056 |
| Vishrut Patel | AU2340216 |

## Goal of M3
M2 built the pipeline. M3 runs the full experimental plan under identical
conditions for every method and reports what the evidence shows: our
latent-variable PGM (probabilistic CCA, learned by CCA-MLE and by EM) against
the base paper's baselines (PCA, SCCA, GN-SCCA with HD / PD / OPD deflation),
with uncertainty across seeds and folds, ablations, robustness and error
analysis, and a revised SOTA position.

## Contents
| Folder | Contents |
|---|---|
| `Report/` | `M3_G02_Report.pdf` |
| `Results/` | Result files, KPI tables and plots |

## Data
Unchanged from M2: the processed TCGA-BRCA cohort released with the base paper,
already in this repository under `03_brca/data/splits_journal/`. See
[`M2_G02/Data/README.md`](../M2_G02/Data/README.md).

## Base paper
V. Subramanian, T. Syeda-Mahmood, M. N. Do, "Modelling-based joint embedding of
histology and genomics using canonical correlation analysis for breast cancer
survival prediction," *Artificial Intelligence in Medicine* 149:102787, 2024.
