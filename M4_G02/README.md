# M4_G02: Final Submission

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

## Project
We predict breast cancer survival from two data sources for the same patients:
gene expression (genomics, X) and features extracted from tissue images
(histology, Y). We model them with a latent-variable PGM, probabilistic CCA:

```
z ~ N(0, I)      x | z ~ N(Wx z + mu_x, Psi_x)      y | z ~ N(Wy z + mu_y, Psi_y)
p(z, x, y) = p(z) p(x | z) p(y | z)        so  X ⊥ Y | Z
```

- **Learning:** closed-form maximum likelihood (CCA-MLE) and EM
- **Inference:** exact Gaussian posteriors p(z | x, y), p(z | x), p(z | y); a missing modality is marginalised out
- **Prediction:** Cox proportional hazards on the posterior mean of z, evaluated with the C-index
- **Baselines:** the base paper's PCA, SCCA and GN-SCCA (HD, PD, OPD deflation), rerun from the authors' trained models, and a non-PGM PCA baseline

**Base paper:** V. Subramanian, T. Syeda-Mahmood, M. N. Do, "Modelling-based joint embedding of histology and genomics using canonical correlation analysis for breast cancer survival prediction," *Artificial Intelligence in Medicine* 149:102787, 2024.

## Contents
| Folder | Contents |
|---|---|
| `Report/` | `M4_G02_Report.pdf` |
| `Presentation/` | `M4_G02_Slides.pdf` |
| `Video/` | `M4_G02_Video_Link.txt` |
| `Code/` | Complete source code (model, learning, inference, baseline, experiments, statistics) and its README |
| `Data/` | Dataset source and preparation |
| `Results/` | Final KPI tables, statistical tables, plots and raw result files |

## Results folder
| Path | Contents |
|---|---|
| `Results/SUMMARY.md`, `Results/kpi_table.csv` | KPI summary: baseline reproduction, C-index, likelihood, posterior entropy, latent-dimension sweeps |
| `Results/fig1` to `fig6` | Graphical model, EM convergence, C-index comparisons, latent-dimension sweeps, posterior entropy, comparison with the paper |
| `Results/stats/` | Mean, SD, 95% CI, paired t-test and Wilcoxon |
| `Results/probability_quality/` | NLL, Brier score, integrated Brier score, calibration plot and ECE |
| `Results/robustness/` | Missing modality, observation noise and training-size results |
| `Results/analysis/`, `Results/plots/` | Aggregated tables and plots for all experiments and ablations |
| `Results/*.json` | Raw per-fold results for every run |

## Reproduce
From the repository root:

```bash
pip install -r M4_G02/Code/requirements.txt -r M4_G02/Code/requirements-paper.txt pandas tabulate
python M4_G02/Code/run_all.py .
```

This regenerates everything in `Results/` (about 6 to 10 minutes). Details,
seeds and individual commands are in [`Code/README.md`](Code/README.md); the
dataset is described in [`Data/README.md`](Data/README.md).
