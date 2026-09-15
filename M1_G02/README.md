# M1_G02: Project Proposal and Scoping

**Project:** Breast Cancer Survival Prediction using a Latent-Variable PGM
**Course:** CSE516 Probabilistic Graphical Models, Ahmedabad University
**Track:** Applied AI (Healthcare)

## Team
| Name | Enrollment No. |
|---|---|
| Dhanya Achari | AU2420163 |
| Shifataj Shaikh | AU2420174 |
| Paarth Jawaharani | AU2420056 |
| Vishrut Patel | AU2340216 |

## Summary
We model genomic (X) and histology-derived (Y) data from TCGA-BRCA through a
shared latent variable Z (Z → X, Z → Y). The CCA-based parameter learning of
Subramanian et al. (2024) is reproduced as the baseline, and EM-based parameter
learning is proposed for the same PGM. We perform latent-state inference with
P(Z|X,Y), P(Z|X) and P(Z|Y), and evaluate survival prediction with the C-index.

## Folder Contents
| Folder | Contents |
|---|---|
| `Report/` | `M1_G02_Report.pdf`: proposal, PGM mapping, SOTA position, baseline, KPIs, plan |
| `Video/` | Link to the M1 video walkthrough |
| `Data/` | Dataset source and description |
| `Code/` | No new code at M1; baseline code is in the repository root |
| `Results/` | No results at M1 |

## Reference Work
V. Subramanian, T. Syeda-Mahmood, and M. N. Do, "Modelling-based joint embedding
of histology and genomics using canonical correlation analysis for breast cancer
survival prediction," *Artificial Intelligence in Medicine*, vol. 149, 102787, 2024.
doi: [10.1016/j.artmed.2024.102787](https://doi.org/10.1016/j.artmed.2024.102787)
