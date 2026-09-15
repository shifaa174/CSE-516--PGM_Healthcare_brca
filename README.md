# CSE516 Probabilistic Graphical Models: Breast Cancer Survival Prediction (Group G02)

**Course:** CSE516 Probabilistic Graphical Models, Ahmedabad University
**Instructor:** Prof. Dhaval Patel
**Track:** Applied AI (Healthcare)

## Team
| Name | Enrollment No. |
|---|---|
| Dhanya Achari | AU2420163 |
| Shifataj Shaikh | AU2420174 |
| Paarth Jawaharani | AU2420056 |
| Vishrut Patel | AU2340216 |

## Project Overview
We study multimodal breast cancer survival prediction on TCGA-BRCA using a
latent-variable PGM (Z → X, Z → Y), where X is genomic data and Y is
histology-derived data. We first reproduce the CCA-based parameter learning of
the reference work as our baseline. We then investigate EM-based parameter
learning for the same model, explicit latent-state inference
(P(Z|X,Y), P(Z|X), P(Z|Y)), and missing-modality settings.
Survival prediction is evaluated with the C-index.

## Milestone Submissions
| Milestone | Folder | Status |
|---|---|---|
| M1: Proposal and Scoping | [`M1_G02/`](M1_G02/) | Submitted |
| M2 | — | Upcoming |
| M3 | — | Upcoming |
| M4 | — | Upcoming |

## Repository Structure
- `M1_G02/`: Milestone 1 submission (report, video, data notes)
- `01_methods/`: CCA, sparse/graph-informed CCA, and deflation implementations (reference code)
- `02_simulations/`: simulation experiments (reference code)
- `03_brca/`: TCGA-BRCA survival prediction pipeline and processed data (reference code)

## Attribution
The code in `01_methods/`, `02_simulations/` and `03_brca/` comes from the
original authors' repository and is used as the baseline for this project:

> V. Subramanian, T. Syeda-Mahmood, and M. N. Do, "Modelling-based joint embedding
> of histology and genomics using canonical correlation analysis for breast cancer
> survival prediction," *Artificial Intelligence in Medicine*, vol. 149, 102787, 2024.
> doi: [10.1016/j.artmed.2024.102787](https://doi.org/10.1016/j.artmed.2024.102787)

All extensions (EM-based learning, latent-state inference, missing-modality
experiments) are the work of Group G02.
