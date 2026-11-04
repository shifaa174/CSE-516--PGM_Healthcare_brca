# M3 results aggregation

Generated from existing JSON result files. No experiment was rerun.

## Interpretation cautions
- Each JSON summary reports mean and SD across the five fixed folds.
- Seed runs 42-46 change EM initialisation only. CCA-MLE and baseline results are deterministic repeats.
- `sd_of_run_means` measures variation among JSON run means; it is not a patient-level confidence interval.
- Paired tests compare CCA-MLE and EM across five folds within one JSON run. With only five folds, treat p-values as exploratory, not definitive.
- Do not treat the five seeds as five independent datasets.

## Files
- `run_level_metrics.csv`: one row per result file and metric.
- `fold_level_metrics.csv`: per-fold observations.
- `experiment_summary.csv`: summary across JSON files grouped by experiment tag.
- `cindex_summary_percent.csv`: C-index summary in percentage points.
- `paired_tests_exploratory.csv`: paired fold-level tests.

## C-index summary

| experiment     | metric                         |   n_json_runs |   mean_pct |   sd_between_runs_pct |   mean_fold_sd |
|:---------------|:-------------------------------|--------------:|-----------:|----------------------:|---------------:|
| abl_eminit_cca | c_index/baseline_pca_concat    |             5 |    6608.2  |                0      |         6.966  |
| abl_eminit_cca | c_index/baseline_pca_genomics  |             5 |    6782.4  |                0      |         7.383  |
| abl_eminit_cca | c_index/baseline_pca_histology |             5 |    5763.6  |                0      |         8.49   |
| abl_eminit_cca | c_index/cca-mle/genomics_only  |             5 |    6730.8  |                0      |         8.613  |
| abl_eminit_cca | c_index/cca-mle/histology_only |             5 |    5733.9  |                0      |         6.455  |
| abl_eminit_cca | c_index/cca-mle/multimodal     |             5 |    6543.3  |                0      |         6.581  |
| abl_eminit_cca | c_index/em/genomics_only       |             5 |    6617.7  |                0      |         6.682  |
| abl_eminit_cca | c_index/em/histology_only      |             5 |    5651.8  |                0      |         7.786  |
| abl_eminit_cca | c_index/em/multimodal          |             5 |    6230.7  |                0      |         4.222  |
| abl_noise_diag | c_index/baseline_pca_concat    |             5 |    6608.2  |                0      |         6.966  |
| abl_noise_diag | c_index/baseline_pca_genomics  |             5 |    6782.4  |                0      |         7.383  |
| abl_noise_diag | c_index/baseline_pca_histology |             5 |    5763.6  |                0      |         8.49   |
| abl_noise_diag | c_index/cca-mle/genomics_only  |             5 |    7105    |                0      |         8.605  |
| abl_noise_diag | c_index/cca-mle/histology_only |             5 |    5693.1  |                0      |         8.781  |
| abl_noise_diag | c_index/cca-mle/multimodal     |             5 |    6797.3  |                0      |         6.758  |
| abl_noise_diag | c_index/em/genomics_only       |             5 |    6846.52 |               34.2051 |         8.6556 |
| abl_noise_diag | c_index/em/histology_only      |             5 |    6256.3  |               19.0955 |         7.5416 |
| abl_noise_diag | c_index/em/multimodal          |             5 |    6438.82 |               30.4524 |         2.9438 |
| abl_pen_0.01   | c_index/baseline_pca_concat    |             5 |    6513.5  |                0      |         7.722  |
| abl_pen_0.01   | c_index/baseline_pca_genomics  |             5 |    6703.7  |                0      |         7.436  |
| abl_pen_0.01   | c_index/baseline_pca_histology |             5 |    5832.8  |                0      |         6.696  |
| abl_pen_0.01   | c_index/cca-mle/genomics_only  |             5 |    6702.7  |                0      |         9.684  |
| abl_pen_0.01   | c_index/cca-mle/histology_only |             5 |    5939.3  |                0      |         5.131  |
| abl_pen_0.01   | c_index/cca-mle/multimodal     |             5 |    6549.6  |                0      |         7.139  |
| abl_pen_0.01   | c_index/em/genomics_only       |             5 |    6785.02 |               43.0207 |         8.1908 |
| abl_pen_0.01   | c_index/em/histology_only      |             5 |    5776.98 |               30.2993 |         7.213  |
| abl_pen_0.01   | c_index/em/multimodal          |             5 |    6220.66 |              172.554  |         7.2242 |
| abl_pen_1.0    | c_index/baseline_pca_concat    |             5 |    6591.9  |                0      |         6.047  |
| abl_pen_1.0    | c_index/baseline_pca_genomics  |             5 |    6792.6  |                0      |         6.982  |
| abl_pen_1.0    | c_index/baseline_pca_histology |             5 |    5525.7  |                0      |         8.696  |
| abl_pen_1.0    | c_index/cca-mle/genomics_only  |             5 |    6754    |                0      |         8.047  |
| abl_pen_1.0    | c_index/cca-mle/histology_only |             5 |    5621.7  |                0      |         6.827  |
| abl_pen_1.0    | c_index/cca-mle/multimodal     |             5 |    6435.3  |                0      |         6.637  |
| abl_pen_1.0    | c_index/em/genomics_only       |             5 |    6848.62 |               55.9002 |         8.2438 |
| abl_pen_1.0    | c_index/em/histology_only      |             5 |    5521.74 |               82.5798 |        10.132  |
| abl_pen_1.0    | c_index/em/multimodal          |             5 |    6231.78 |              105.77   |         6.6274 |
| abl_reg_1e-1   | c_index/baseline_pca_concat    |             5 |    6608.2  |                0      |         6.966  |
| abl_reg_1e-1   | c_index/baseline_pca_genomics  |             5 |    6782.4  |                0      |         7.383  |
| abl_reg_1e-1   | c_index/baseline_pca_histology |             5 |    5763.6  |                0      |         8.49   |
| abl_reg_1e-1   | c_index/cca-mle/genomics_only  |             5 |    6823.8  |                0      |         9.09   |
| abl_reg_1e-1   | c_index/cca-mle/histology_only |             5 |    5801.3  |                0      |         6.392  |
| abl_reg_1e-1   | c_index/cca-mle/multimodal     |             5 |    6577.9  |                0      |         7.042  |
| abl_reg_1e-1   | c_index/em/genomics_only       |             5 |    6912.98 |               50.5729 |         8.155  |
| abl_reg_1e-1   | c_index/em/histology_only      |             5 |    5673.82 |               57.0304 |         8.5552 |
| abl_reg_1e-1   | c_index/em/multimodal          |             5 |    6333.5  |              128.556  |         6.8396 |
| abl_reg_1e-3   | c_index/baseline_pca_concat    |             5 |    6608.2  |                0      |         6.966  |
| abl_reg_1e-3   | c_index/baseline_pca_genomics  |             5 |    6782.4  |                0      |         7.383  |
| abl_reg_1e-3   | c_index/baseline_pca_histology |             5 |    5763.6  |                0      |         8.49   |
| abl_reg_1e-3   | c_index/cca-mle/genomics_only  |             5 |    6644.8  |                0      |         8.452  |
| abl_reg_1e-3   | c_index/cca-mle/histology_only |             5 |    5777    |                0      |         7.128  |
| abl_reg_1e-3   | c_index/cca-mle/multimodal     |             5 |    6434.2  |                0      |         6.614  |
| abl_reg_1e-3   | c_index/em/genomics_only       |             5 |    6912.98 |               50.5729 |         8.155  |
| abl_reg_1e-3   | c_index/em/histology_only      |             5 |    5673.82 |               57.0304 |         8.5552 |
| abl_reg_1e-3   | c_index/em/multimodal          |             5 |    6333.5  |              128.556  |         6.8396 |
| paperpca       | c_index/baseline_pca_concat    |             5 |    6608.2  |                0      |         6.966  |
| paperpca       | c_index/baseline_pca_genomics  |             5 |    6782.4  |                0      |         7.383  |
| paperpca       | c_index/baseline_pca_histology |             5 |    5763.6  |                0      |         8.49   |
| paperpca       | c_index/cca-mle/genomics_only  |             5 |    6730.8  |                0      |         8.613  |
| paperpca       | c_index/cca-mle/histology_only |             5 |    5733.9  |                0      |         6.455  |
| paperpca       | c_index/cca-mle/multimodal     |             5 |    6543.3  |                0      |         6.581  |
| paperpca       | c_index/em/genomics_only       |             5 |    6912.98 |               50.5729 |         8.155  |
| paperpca       | c_index/em/histology_only      |             5 |    5673.82 |               57.0304 |         8.5552 |
| paperpca       | c_index/em/multimodal          |             5 |    6333.5  |              128.556  |         6.8396 |

## Paired tests (exploratory)

| file                                  | experiment     |   seed | outcome        |   n_paired_folds |   cca_mean |   em_mean |   mean_em_minus_cca |   wilcoxon_p_exploratory |   paired_t_p_exploratory |
|:--------------------------------------|:---------------|-------:|:---------------|-----------------:|-----------:|----------:|--------------------:|-------------------------:|-------------------------:|
| m2_results_abl_eminit_cca_seed42.json | abl_eminit_cca |     42 | multimodal     |                5 |   0.654335 |  0.623065 |        -0.0312699   |                   0.375  |                0.374981  |
| m2_results_abl_eminit_cca_seed42.json | abl_eminit_cca |     42 | genomics_only  |                5 |   0.67308  |  0.66177  |        -0.0113095   |                   1      |                0.765737  |
| m2_results_abl_eminit_cca_seed42.json | abl_eminit_cca |     42 | histology_only |                5 |   0.573392 |  0.565179 |        -0.0082127   |                   0.8125 |                0.74911   |
| m2_results_abl_eminit_cca_seed43.json | abl_eminit_cca |     43 | multimodal     |                5 |   0.654335 |  0.623065 |        -0.0312699   |                   0.375  |                0.374981  |
| m2_results_abl_eminit_cca_seed43.json | abl_eminit_cca |     43 | genomics_only  |                5 |   0.67308  |  0.66177  |        -0.0113095   |                   1      |                0.765737  |
| m2_results_abl_eminit_cca_seed43.json | abl_eminit_cca |     43 | histology_only |                5 |   0.573392 |  0.565179 |        -0.0082127   |                   0.8125 |                0.74911   |
| m2_results_abl_eminit_cca_seed44.json | abl_eminit_cca |     44 | multimodal     |                5 |   0.654335 |  0.623065 |        -0.0312699   |                   0.375  |                0.374981  |
| m2_results_abl_eminit_cca_seed44.json | abl_eminit_cca |     44 | genomics_only  |                5 |   0.67308  |  0.66177  |        -0.0113095   |                   1      |                0.765737  |
| m2_results_abl_eminit_cca_seed44.json | abl_eminit_cca |     44 | histology_only |                5 |   0.573392 |  0.565179 |        -0.0082127   |                   0.8125 |                0.74911   |
| m2_results_abl_eminit_cca_seed45.json | abl_eminit_cca |     45 | multimodal     |                5 |   0.654335 |  0.623065 |        -0.0312699   |                   0.375  |                0.374981  |
| m2_results_abl_eminit_cca_seed45.json | abl_eminit_cca |     45 | genomics_only  |                5 |   0.67308  |  0.66177  |        -0.0113095   |                   1      |                0.765737  |
| m2_results_abl_eminit_cca_seed45.json | abl_eminit_cca |     45 | histology_only |                5 |   0.573392 |  0.565179 |        -0.0082127   |                   0.8125 |                0.74911   |
| m2_results_abl_eminit_cca_seed46.json | abl_eminit_cca |     46 | multimodal     |                5 |   0.654335 |  0.623065 |        -0.0312699   |                   0.375  |                0.374981  |
| m2_results_abl_eminit_cca_seed46.json | abl_eminit_cca |     46 | genomics_only  |                5 |   0.67308  |  0.66177  |        -0.0113095   |                   1      |                0.765737  |
| m2_results_abl_eminit_cca_seed46.json | abl_eminit_cca |     46 | histology_only |                5 |   0.573392 |  0.565179 |        -0.0082127   |                   0.8125 |                0.74911   |
| m2_results_abl_noise_diag_seed42.json | abl_noise_diag |     42 | multimodal     |                5 |   0.679727 |  0.640212 |        -0.0395149   |                   0.3125 |                0.231613  |
| m2_results_abl_noise_diag_seed42.json | abl_noise_diag |     42 | genomics_only  |                5 |   0.710495 |  0.684321 |        -0.0261748   |                   0.3125 |                0.31165   |
| m2_results_abl_noise_diag_seed42.json | abl_noise_diag |     42 | histology_only |                5 |   0.569315 |  0.62334  |         0.0540254   |                   0.1875 |                0.129434  |
| m2_results_abl_noise_diag_seed43.json | abl_noise_diag |     43 | multimodal     |                5 |   0.679727 |  0.647262 |        -0.0324642   |                   0.4375 |                0.352136  |
| m2_results_abl_noise_diag_seed43.json | abl_noise_diag |     43 | genomics_only  |                5 |   0.710495 |  0.681004 |        -0.0294912   |                   0.3125 |                0.239166  |
| m2_results_abl_noise_diag_seed43.json | abl_noise_diag |     43 | histology_only |                5 |   0.569315 |  0.627384 |         0.0580696   |                   0.125  |                0.0701006 |
| m2_results_abl_noise_diag_seed44.json | abl_noise_diag |     44 | multimodal     |                5 |   0.679727 |  0.644786 |        -0.0349407   |                   0.3125 |                0.326129  |
| m2_results_abl_noise_diag_seed44.json | abl_noise_diag |     44 | genomics_only  |                5 |   0.710495 |  0.688835 |        -0.0216603   |                   0.625  |                0.352235  |
| m2_results_abl_noise_diag_seed44.json | abl_noise_diag |     44 | histology_only |                5 |   0.569315 |  0.62478  |         0.0554655   |                   0.125  |                0.0867292 |
| m2_results_abl_noise_diag_seed45.json | abl_noise_diag |     45 | multimodal     |                5 |   0.679727 |  0.641209 |        -0.0385174   |                   0.3125 |                0.258845  |
| m2_results_abl_noise_diag_seed45.json | abl_noise_diag |     45 | genomics_only  |                5 |   0.710495 |  0.681741 |        -0.028754    |                   0.3125 |                0.252305  |
| m2_results_abl_noise_diag_seed45.json | abl_noise_diag |     45 | histology_only |                5 |   0.569315 |  0.627835 |         0.0585209   |                   0.125  |                0.101003  |
| m2_results_abl_noise_diag_seed46.json | abl_noise_diag |     46 | multimodal     |                5 |   0.679727 |  0.645945 |        -0.0337821   |                   0.4375 |                0.325337  |
| m2_results_abl_noise_diag_seed46.json | abl_noise_diag |     46 | genomics_only  |                5 |   0.710495 |  0.687361 |        -0.0231344   |                   0.3125 |                0.370883  |
| m2_results_abl_noise_diag_seed46.json | abl_noise_diag |     46 | histology_only |                5 |   0.569315 |  0.624806 |         0.0554919   |                   0.1875 |                0.125525  |
| m2_results_abl_pen_0.01_seed42.json   | abl_pen_0.01   |     42 | multimodal     |                5 |   0.654958 |  0.621799 |        -0.0331586   |                   0.3125 |                0.217147  |
| m2_results_abl_pen_0.01_seed42.json   | abl_pen_0.01   |     42 | genomics_only  |                5 |   0.670273 |  0.685344 |         0.0150708   |                   0.4375 |                0.463507  |
| m2_results_abl_pen_0.01_seed42.json   | abl_pen_0.01   |     42 | histology_only |                5 |   0.593932 |  0.576055 |        -0.0178775   |                   0.4375 |                0.535582  |
| m2_results_abl_pen_0.01_seed43.json   | abl_pen_0.01   |     43 | multimodal     |                5 |   0.654958 |  0.633694 |        -0.0212637   |                   0.4375 |                0.411248  |
| m2_results_abl_pen_0.01_seed43.json   | abl_pen_0.01   |     43 | genomics_only  |                5 |   0.670273 |  0.675199 |         0.00492509  |                   0.8125 |                0.838249  |
| m2_results_abl_pen_0.01_seed43.json   | abl_pen_0.01   |     43 | histology_only |                5 |   0.593932 |  0.575549 |        -0.018383    |                   0.4375 |                0.520554  |
| m2_results_abl_pen_0.01_seed44.json   | abl_pen_0.01   |     44 | multimodal     |                5 |   0.654958 |  0.628655 |        -0.0263025   |                   0.1875 |                0.171216  |
| m2_results_abl_pen_0.01_seed44.json   | abl_pen_0.01   |     44 | genomics_only  |                5 |   0.670273 |  0.679794 |         0.00952048  |                   0.625  |                0.650324  |
| m2_results_abl_pen_0.01_seed44.json   | abl_pen_0.01   |     44 | histology_only |                5 |   0.593932 |  0.580144 |        -0.0137884   |                   0.625  |                0.591294  |
| m2_results_abl_pen_0.01_seed45.json   | abl_pen_0.01   |     45 | multimodal     |                5 |   0.654958 |  0.592464 |        -0.0624934   |                   0.125  |                0.0655918 |
| m2_results_abl_pen_0.01_seed45.json   | abl_pen_0.01   |     45 | genomics_only  |                5 |   0.670273 |  0.677335 |         0.00706127  |                   1      |                0.760983  |
| m2_results_abl_pen_0.01_seed45.json   | abl_pen_0.01   |     45 | histology_only |                5 |   0.593932 |  0.575023 |        -0.0189096   |                   0.625  |                0.528188  |
| m2_results_abl_pen_0.01_seed46.json   | abl_pen_0.01   |     46 | multimodal     |                5 |   0.654958 |  0.633716 |        -0.0212414   |                   0.4375 |                0.387049  |
| m2_results_abl_pen_0.01_seed46.json   | abl_pen_0.01   |     46 | genomics_only  |                5 |   0.670273 |  0.674849 |         0.00457505  |                   0.8125 |                0.838201  |
| m2_results_abl_pen_0.01_seed46.json   | abl_pen_0.01   |     46 | histology_only |                5 |   0.593932 |  0.581731 |        -0.0122016   |                   0.8125 |                0.613039  |
| m2_results_abl_pen_1.0_seed42.json    | abl_pen_1.0    |     42 | multimodal     |                5 |   0.643527 |  0.619707 |        -0.0238198   |                   0.4375 |                0.362833  |
| m2_results_abl_pen_1.0_seed42.json    | abl_pen_1.0    |     42 | genomics_only  |                5 |   0.675397 |  0.689389 |         0.0139915   |                   0.8125 |                0.38758   |
| m2_results_abl_pen_1.0_seed42.json    | abl_pen_1.0    |     42 | histology_only |                5 |   0.562171 |  0.553866 |        -0.00830502  |                   0.625  |                0.722219  |
| m2_results_abl_pen_1.0_seed43.json    | abl_pen_1.0    |     43 | multimodal     |                5 |   0.643527 |  0.625037 |        -0.0184902   |                   0.4375 |                0.481957  |
| m2_results_abl_pen_1.0_seed43.json    | abl_pen_1.0    |     43 | genomics_only  |                5 |   0.675397 |  0.684204 |         0.00880728  |                   0.625  |                0.552816  |
| m2_results_abl_pen_1.0_seed43.json    | abl_pen_1.0    |     43 | histology_only |                5 |   0.562171 |  0.553952 |        -0.00821848  |                   0.625  |                0.749786  |
| m2_results_abl_pen_1.0_seed44.json    | abl_pen_1.0    |     44 | multimodal     |                5 |   0.643527 |  0.628482 |        -0.015045    |                   0.4375 |                0.534251  |
| m2_results_abl_pen_1.0_seed44.json    | abl_pen_1.0    |     44 | genomics_only  |                5 |   0.675397 |  0.676093 |         0.000696055 |                   1      |                0.956286  |
| m2_results_abl_pen_1.0_seed44.json    | abl_pen_1.0    |     44 | histology_only |                5 |   0.562171 |  0.553225 |        -0.00894535  |                   0.625  |                0.677438  |
| m2_results_abl_pen_1.0_seed45.json    | abl_pen_1.0    |     45 | multimodal     |                5 |   0.643527 |  0.607238 |        -0.0362897   |                   0.1875 |                0.158191  |
| m2_results_abl_pen_1.0_seed45.json    | abl_pen_1.0    |     45 | genomics_only  |                5 |   0.675397 |  0.684567 |         0.00917027  |                   0.625  |                0.427954  |
| m2_results_abl_pen_1.0_seed45.json    | abl_pen_1.0    |     45 | histology_only |                5 |   0.562171 |  0.538609 |        -0.0235618   |                   0.4375 |                0.355745  |
| m2_results_abl_pen_1.0_seed46.json    | abl_pen_1.0    |     46 | multimodal     |                5 |   0.643527 |  0.635421 |        -0.00810677  |                   0.8125 |                0.728796  |
| m2_results_abl_pen_1.0_seed46.json    | abl_pen_1.0    |     46 | genomics_only  |                5 |   0.675397 |  0.690058 |         0.0146613   |                   0.3125 |                0.320025  |
| m2_results_abl_pen_1.0_seed46.json    | abl_pen_1.0    |     46 | histology_only |                5 |   0.562171 |  0.561207 |        -0.000963672 |                   1      |                0.966335  |
| m2_results_abl_reg_1e-1_seed42.json   | abl_reg_1e-1   |     42 | multimodal     |                5 |   0.657786 |  0.634418 |        -0.0233684   |                   0.4375 |                0.344303  |
| m2_results_abl_reg_1e-1_seed42.json   | abl_reg_1e-1   |     42 | genomics_only  |                5 |   0.68238  |  0.699504 |         0.0171237   |                   0.4375 |                0.375929  |
| m2_results_abl_reg_1e-1_seed42.json   | abl_reg_1e-1   |     42 | histology_only |                5 |   0.580129 |  0.569276 |        -0.0108532   |                   0.8125 |                0.662096  |
| m2_results_abl_reg_1e-1_seed43.json   | abl_reg_1e-1   |     43 | multimodal     |                5 |   0.657786 |  0.635763 |        -0.0220231   |                   0.4375 |                0.441612  |
| m2_results_abl_reg_1e-1_seed43.json   | abl_reg_1e-1   |     43 | genomics_only  |                5 |   0.68238  |  0.688719 |         0.00633886  |                   0.8125 |                0.743156  |
| m2_results_abl_reg_1e-1_seed43.json   | abl_reg_1e-1   |     43 | histology_only |                5 |   0.580129 |  0.565287 |        -0.0148423   |                   0.625  |                0.592696  |
| m2_results_abl_reg_1e-1_seed44.json   | abl_reg_1e-1   |     44 | multimodal     |                5 |   0.657786 |  0.641112 |        -0.0166745   |                   0.625  |                0.47354   |
| m2_results_abl_reg_1e-1_seed44.json   | abl_reg_1e-1   |     44 | genomics_only  |                5 |   0.68238  |  0.686249 |         0.00386893  |                   1      |                0.825593  |
| m2_results_abl_reg_1e-1_seed44.json   | abl_reg_1e-1   |     44 | histology_only |                5 |   0.580129 |  0.566549 |        -0.0135804   |                   0.625  |                0.56819   |
| m2_results_abl_reg_1e-1_seed45.json   | abl_reg_1e-1   |     45 | multimodal     |                5 |   0.657786 |  0.611443 |        -0.0463432   |                   0.1875 |                0.169432  |
| m2_results_abl_reg_1e-1_seed45.json   | abl_reg_1e-1   |     45 | genomics_only  |                5 |   0.68238  |  0.692184 |         0.00980441  |                   0.4375 |                0.553473  |
| m2_results_abl_reg_1e-1_seed45.json   | abl_reg_1e-1   |     45 | histology_only |                5 |   0.580129 |  0.56012  |        -0.0200088   |                   0.625  |                0.454838  |
| m2_results_abl_reg_1e-1_seed46.json   | abl_reg_1e-1   |     46 | multimodal     |                5 |   0.657786 |  0.644017 |        -0.0137691   |                   0.625  |                0.533964  |
| m2_results_abl_reg_1e-1_seed46.json   | abl_reg_1e-1   |     46 | genomics_only  |                5 |   0.68238  |  0.689838 |         0.0074575   |                   0.8125 |                0.648735  |
| m2_results_abl_reg_1e-1_seed46.json   | abl_reg_1e-1   |     46 | histology_only |                5 |   0.580129 |  0.575665 |        -0.00446398  |                   1      |                0.839219  |
| m2_results_abl_reg_1e-3_seed42.json   | abl_reg_1e-3   |     42 | multimodal     |                5 |   0.643423 |  0.634418 |        -0.00900528  |                   0.625  |                0.687181  |
| m2_results_abl_reg_1e-3_seed42.json   | abl_reg_1e-3   |     42 | genomics_only  |                5 |   0.664483 |  0.699504 |         0.0350203   |                   0.0625 |                0.0421344 |
| m2_results_abl_reg_1e-3_seed42.json   | abl_reg_1e-3   |     42 | histology_only |                5 |   0.577701 |  0.569276 |        -0.00842491  |                   1      |                0.764615  |
| m2_results_abl_reg_1e-3_seed43.json   | abl_reg_1e-3   |     43 | multimodal     |                5 |   0.643423 |  0.635763 |        -0.00765998  |                   0.8125 |                0.77541   |
| m2_results_abl_reg_1e-3_seed43.json   | abl_reg_1e-3   |     43 | genomics_only  |                5 |   0.664483 |  0.688719 |         0.0242355   |                   0.1875 |                0.132591  |
| m2_results_abl_reg_1e-3_seed43.json   | abl_reg_1e-3   |     43 | histology_only |                5 |   0.577701 |  0.565287 |        -0.012414    |                   0.625  |                0.661633  |
| m2_results_abl_reg_1e-3_seed44.json   | abl_reg_1e-3   |     44 | multimodal     |                5 |   0.643423 |  0.641112 |        -0.00231132  |                   1      |                0.912104  |
| m2_results_abl_reg_1e-3_seed44.json   | abl_reg_1e-3   |     44 | genomics_only  |                5 |   0.664483 |  0.686249 |         0.0217656   |                   0.3125 |                0.136581  |
| m2_results_abl_reg_1e-3_seed44.json   | abl_reg_1e-3   |     44 | histology_only |                5 |   0.577701 |  0.566549 |        -0.0111521   |                   1      |                0.644812  |
| m2_results_abl_reg_1e-3_seed45.json   | abl_reg_1e-3   |     45 | multimodal     |                5 |   0.643423 |  0.611443 |        -0.03198     |                   0.3125 |                0.254574  |
| m2_results_abl_reg_1e-3_seed45.json   | abl_reg_1e-3   |     45 | genomics_only  |                5 |   0.664483 |  0.692184 |         0.0277011   |                   0.125  |                0.176194  |
| m2_results_abl_reg_1e-3_seed45.json   | abl_reg_1e-3   |     45 | histology_only |                5 |   0.577701 |  0.56012  |        -0.0175805   |                   0.4375 |                0.476538  |
| m2_results_abl_reg_1e-3_seed46.json   | abl_reg_1e-3   |     46 | multimodal     |                5 |   0.643423 |  0.644017 |         0.00059408  |                   1      |                0.978074  |
| m2_results_abl_reg_1e-3_seed46.json   | abl_reg_1e-3   |     46 | genomics_only  |                5 |   0.664483 |  0.689838 |         0.0253542   |                   0.0625 |                0.0336397 |
| m2_results_abl_reg_1e-3_seed46.json   | abl_reg_1e-3   |     46 | histology_only |                5 |   0.577701 |  0.575665 |        -0.00203568  |                   1      |                0.937937  |
| m2_results_paperpca_seed42.json       | paperpca       |     42 | multimodal     |                5 |   0.654335 |  0.634418 |        -0.0199173   |                   0.4375 |                0.468873  |
| m2_results_paperpca_seed42.json       | paperpca       |     42 | genomics_only  |                5 |   0.67308  |  0.699504 |         0.026424    |                   0.3125 |                0.179877  |
| m2_results_paperpca_seed42.json       | paperpca       |     42 | histology_only |                5 |   0.573392 |  0.569276 |        -0.0041158   |                   1      |                0.865628  |
| m2_results_paperpca_seed43.json       | paperpca       |     43 | multimodal     |                5 |   0.654335 |  0.635763 |        -0.018572    |                   0.4375 |                0.557852  |
| m2_results_paperpca_seed43.json       | paperpca       |     43 | genomics_only  |                5 |   0.67308  |  0.688719 |         0.0156392   |                   0.8125 |                0.362008  |
| m2_results_paperpca_seed43.json       | paperpca       |     43 | histology_only |                5 |   0.573392 |  0.565287 |        -0.00810489  |                   0.8125 |                0.750191  |
| m2_results_paperpca_seed44.json       | paperpca       |     44 | multimodal     |                5 |   0.654335 |  0.641112 |        -0.0132234   |                   0.625  |                0.602306  |
| m2_results_paperpca_seed44.json       | paperpca       |     44 | genomics_only  |                5 |   0.67308  |  0.686249 |         0.0131693   |                   0.625  |                0.416238  |
| m2_results_paperpca_seed44.json       | paperpca       |     44 | histology_only |                5 |   0.573392 |  0.566549 |        -0.00684295  |                   1      |                0.743392  |
| m2_results_paperpca_seed45.json       | paperpca       |     45 | multimodal     |                5 |   0.654335 |  0.611443 |        -0.0428921   |                   0.1875 |                0.193273  |
| m2_results_paperpca_seed45.json       | paperpca       |     45 | genomics_only  |                5 |   0.67308  |  0.692184 |         0.0191048   |                   0.3125 |                0.282328  |
| m2_results_paperpca_seed45.json       | paperpca       |     45 | histology_only |                5 |   0.573392 |  0.56012  |        -0.0132714   |                   0.625  |                0.550083  |
| m2_results_paperpca_seed46.json       | paperpca       |     46 | multimodal     |                5 |   0.654335 |  0.644017 |        -0.010318    |                   0.625  |                0.679927  |
| m2_results_paperpca_seed46.json       | paperpca       |     46 | genomics_only  |                5 |   0.67308  |  0.689838 |         0.0167579   |                   0.3125 |                0.28057   |
| m2_results_paperpca_seed46.json       | paperpca       |     46 | histology_only |                5 |   0.573392 |  0.575665 |         0.00227343  |                   1      |                0.918516  |