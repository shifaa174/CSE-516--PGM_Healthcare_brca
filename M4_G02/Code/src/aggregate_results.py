import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

METRICS = [
    "c_index/cca-mle/multimodal", "c_index/em/multimodal",
    "c_index/baseline_pca_concat",
    "c_index/cca-mle/genomics_only", "c_index/em/genomics_only",
    "c_index/baseline_pca_genomics",
    "c_index/cca-mle/histology_only", "c_index/em/histology_only",
    "c_index/baseline_pca_histology",
    "loglik_per_sample/cca-mle/test", "loglik_per_sample/em/test",
    "fit_seconds/cca-mle", "fit_seconds/em",
    "n_iter/cca-mle", "n_iter/em",
]

def experiment_group(path, cfg):
    tag = str(cfg.get("tag") or path.stem)
    return re.sub(r"_seed\d+$", "", tag)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="M4_G02/Results")
    ap.add_argument("--out", default="M4_G02/Results/analysis")
    args = ap.parse_args()
    src, out = Path(args.results), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    run_rows, fold_rows, paired_rows = [], [], []
    for path in sorted(src.glob("*.json")):
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(f"SKIP {path.name}: {exc}")
            continue
        cfg = obj.get("config", {})
        if not obj.get("per_fold") or not obj.get("summary"):
            continue
        group = experiment_group(path, cfg)
        if not re.search(r"_seed\d+$", str(cfg.get("tag") or path.stem)):
            continue
        seed = cfg.get("seed", np.nan)
        for metric in METRICS:
            summary = obj["summary"].get(metric)
            if summary:
                run_rows.append({
                    "file": path.name, "experiment": group, "seed": seed,
                    "variant": cfg.get("variant", ""),
                    "noise": cfg.get("noise", ""),
                    "em_init": cfg.get("em_init", ""),
                    "reg": cfg.get("reg", np.nan),
                    "cox_penalizer": cfg.get("penalizer", np.nan),
                    "metric": metric, "mean_over_folds": summary.get("mean"),
                    "sd_over_folds": summary.get("std"),
                })
        for fold in obj["per_fold"]:
            base = {"file": path.name, "experiment": group, "seed": seed,
                    "fold": fold.get("fold")}
            for metric in METRICS:
                parts = metric.split("/")
                value = None
                if parts[0] == "c_index" and len(parts) == 3:
                    if parts[1] in ("cca-mle", "em"):
                        value = fold.get(parts[1], {}).get(parts[2], {}).get("c_index")
                elif parts[0] == "c_index" and len(parts) == 2:
                    value = fold.get(parts[1])
                elif metric.startswith("loglik_per_sample/"):
                    method, split = parts[1], parts[2]
                    value = fold.get(method, {}).get(split + "_loglik_per_sample")
                elif metric.startswith("fit_seconds/"):
                    value = fold.get(parts[1], {}).get("fit_seconds")
                elif metric.startswith("n_iter/"):
                    value = fold.get(parts[1], {}).get("n_iter")
                if value is not None:
                    fold_rows.append({**base, "metric": metric, "value": value})

        for outcome in ("multimodal", "genomics_only", "histology_only"):
            a = [f.get("cca-mle", {}).get(outcome, {}).get("c_index") for f in obj["per_fold"]]
            b = [f.get("em", {}).get(outcome, {}).get("c_index") for f in obj["per_fold"]]
            pairs = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
            if len(pairs) >= 3:
                x, y = np.array(pairs).T
                try:
                    wil = stats.wilcoxon(x, y, zero_method="wilcox", alternative="two-sided")
                    wp = float(wil.pvalue)
                except Exception:
                    wp = np.nan
                try:
                    tp = float(stats.ttest_rel(x, y).pvalue)
                except Exception:
                    tp = np.nan
                paired_rows.append({
                    "file": path.name, "experiment": group, "seed": seed,
                    "outcome": outcome, "n_paired_folds": len(pairs),
                    "cca_mean": float(x.mean()), "em_mean": float(y.mean()),
                    "mean_em_minus_cca": float((y-x).mean()),
                    "wilcoxon_p_exploratory": wp, "paired_t_p_exploratory": tp,
                })

    runs = pd.DataFrame(run_rows)
    folds = pd.DataFrame(fold_rows)
    paired = pd.DataFrame(paired_rows)
    if runs.empty:
        raise SystemExit(f"No compatible result JSONs found in {src}")

    runs.to_csv(out / "run_level_metrics.csv", index=False)
    folds.to_csv(out / "fold_level_metrics.csv", index=False)
    paired.to_csv(out / "paired_tests_exploratory.csv", index=False)

    agg = (runs.groupby(["experiment", "metric"], dropna=False)
           .agg(n_json_runs=("file", "nunique"),
                mean_of_run_means=("mean_over_folds", "mean"),
                sd_of_run_means=("mean_over_folds", "std"),
                mean_fold_sd=("sd_over_folds", "mean"))
           .reset_index())
    agg.to_csv(out / "experiment_summary.csv", index=False)

    c = agg[agg.metric.str.startswith("c_index/")].copy()
    c["mean_pct"] = c["mean_of_run_means"]
    c["sd_between_runs_pct"] = c["sd_of_run_means"]
    c.to_csv(out / "cindex_summary_percent.csv", index=False)

    readme = [
        "# Results aggregation",
        "",
        "Generated from existing JSON result files. No experiment was rerun.",
        "",
        "## Interpretation cautions",
        "- Each JSON summary reports mean and SD across the five fixed folds.",
        "- Seed runs 42-46 change EM initialisation only. CCA-MLE and baseline results are deterministic repeats.",
        "- `sd_of_run_means` measures variation among JSON run means; it is not a patient-level confidence interval.",
        "- Paired tests compare CCA-MLE and EM across five folds within one JSON run. With only five folds, treat p-values as exploratory, not definitive.",
        "- Do not treat the five seeds as five independent datasets.",
        "",
        "## Files",
        "- `run_level_metrics.csv`: one row per result file and metric.",
        "- `fold_level_metrics.csv`: per-fold observations.",
        "- `experiment_summary.csv`: summary across JSON files grouped by experiment tag.",
        "- `cindex_summary_percent.csv`: C-index summary in percentage points.",
        "- `paired_tests_exploratory.csv`: paired fold-level tests.",
        "",
        "## C-index summary",
        "",
        c[["experiment","metric","n_json_runs","mean_pct","sd_between_runs_pct","mean_fold_sd"]].to_markdown(index=False),
        "",
        "## Paired tests (exploratory)",
        "",
        paired.to_markdown(index=False) if not paired.empty else "No paired comparisons available.",
    ]
    (out / "ANALYSIS_README.md").write_text("\n".join(readme), encoding="utf-8")
    print(f"Saved analysis tables to: {out.resolve()}")

if __name__ == "__main__":
    main()
