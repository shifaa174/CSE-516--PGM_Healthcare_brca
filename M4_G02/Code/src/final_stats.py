import argparse, glob, json, os
import numpy as np
import pandas as pd
from scipy import stats

SETTINGS = ["multimodal", "genomics_only", "histology_only"]
BASELINE = {"multimodal": "baseline_pca_concat", "genomics_only": "baseline_pca_genomics",
            "histology_only": "baseline_pca_histology"}
PAPER_SETTING = {"multimodal": "concat", "genomics_only": "genomics", "histology_only": "imaging"}


def fold_values(results, tag):
    files = sorted(glob.glob(os.path.join(results, f"results_{tag}_seed*.json")))
    if not files:
        return None
    runs = [json.load(open(f, encoding="utf-8"))["per_fold"] for f in files]
    out = {}
    for s in SETTINGS:
        for m in ("cca-mle", "em"):
            out[f"{m}/{s}"] = np.mean([[100 * f[m][s]["c_index"] for f in r] for r in runs], axis=0)
        out[f"baseline/{s}"] = np.mean([[100 * f[BASELINE[s]] for f in r] for r in runs], axis=0)
    for m in ("cca-mle", "em"):
        out[f"{m}/test_nll"] = np.mean([[-f[m]["test_loglik_per_sample"] for f in r] for r in runs], axis=0)
        out[f"{m}/fit_seconds"] = np.mean([[f[m]["fit_seconds"] for f in r] for r in runs], axis=0)
    out["n_seeds"] = len(files)
    return out


def paper_values(results):
    rows = json.load(open(os.path.join(results, "paper_baseline.json"), encoding="utf-8"))["rows"]
    out = {}
    for r in rows:
        if r["rerun"] and r["rerun"].get("per_fold"):
            out[f"paper {r['method']}/{r['setting']}"] = np.array(r["rerun"]["per_fold"], dtype=float)
    return out


def summary(name, v):
    n = len(v)
    m, sd = float(np.mean(v)), float(np.std(v, ddof=1))
    h = stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n)
    return dict(method=name, n_folds=n, mean=round(m, 2), sd=round(sd, 2),
                ci95_low=round(m - h, 2), ci95_high=round(m + h, 2))


def paired(name_a, a, name_b, b):
    d = a - b
    n = len(d)
    sd = float(np.std(d, ddof=1))
    h = stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n)
    t_p = float(stats.ttest_rel(a, b).pvalue)
    w_p = float(stats.wilcoxon(a, b).pvalue) if np.any(d != 0) else 1.0
    return dict(comparison=f"{name_a} - {name_b}", n_folds=n, mean_diff=round(float(d.mean()), 2),
                ci95_low=round(float(d.mean() - h), 2), ci95_high=round(float(d.mean() + h), 2),
                paired_t_p=round(t_p, 4), wilcoxon_p=round(w_p, 4),
                wins_a=int((d > 0).sum()), wins_b=int((d < 0).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="M4_G02/Results")
    ap.add_argument("--out", default="M4_G02/Results/stats")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    block = fold_values(args.results, "paperpca")
    diag = fold_values(args.results, "abl_noise_diag")
    paper = paper_values(args.results)

    rows = []
    for s in SETTINGS:
        rows.append(dict(setting=s, **summary("CCA-MLE (block noise)", block[f"cca-mle/{s}"])))
        rows.append(dict(setting=s, **summary("EM (block noise)", block[f"em/{s}"])))
        if diag:
            rows.append(dict(setting=s, **summary("CCA-MLE (diagonal noise)", diag[f"cca-mle/{s}"])))
            rows.append(dict(setting=s, **summary("EM (diagonal noise)", diag[f"em/{s}"])))
        rows.append(dict(setting=s, **summary("Non-PGM PCA baseline", block[f"baseline/{s}"])))
        for k, v in paper.items():
            if k.endswith("/" + PAPER_SETTING[s]):
                rows.append(dict(setting=s, **summary(k.rsplit("/", 1)[0], v)))
    kpi = pd.DataFrame(rows)
    kpi.to_csv(os.path.join(args.out, "final_cindex_ci.csv"), index=False)

    tests = []
    pcat = paper["paper SCCA/OPD/concat"]
    ppca = paper["paper PCA/concat"]
    for s in SETTINGS:
        cca, em, base = block[f"cca-mle/{s}"], block[f"em/{s}"], block[f"baseline/{s}"]
        tests.append(dict(setting=s, **paired("CCA-MLE", cca, "EM", em)))
        tests.append(dict(setting=s, **paired("CCA-MLE", cca, "Non-PGM PCA baseline", base)))
        tests.append(dict(setting=s, **paired("EM", em, "Non-PGM PCA baseline", base)))
    tests.append(dict(setting="multimodal", **paired("CCA-MLE", block["cca-mle/multimodal"], "paper SCCA/OPD", pcat)))
    tests.append(dict(setting="multimodal", **paired("EM", block["em/multimodal"], "paper SCCA/OPD", pcat)))
    tests.append(dict(setting="multimodal", **paired("CCA-MLE", block["cca-mle/multimodal"], "paper PCA", ppca)))
    if diag:
        tests.append(dict(setting="multimodal", **paired("CCA-MLE diagonal", diag["cca-mle/multimodal"], "CCA-MLE block", block["cca-mle/multimodal"])))
        tests.append(dict(setting="multimodal", **paired("EM diagonal", diag["em/multimodal"], "EM block", block["em/multimodal"])))
        tests.append(dict(setting="multimodal", **paired("CCA-MLE diagonal", diag["cca-mle/multimodal"], "paper SCCA/OPD", pcat)))
    tt = pd.DataFrame(tests)
    tt.to_csv(os.path.join(args.out, "final_paired_tests.csv"), index=False)

    cost = pd.DataFrame([dict(method=m, test_nll_per_sample=round(float(block[f"{k}/test_nll"].mean()), 2),
                              test_nll_sd=round(float(block[f"{k}/test_nll"].std(ddof=1)), 2),
                              fit_seconds=round(float(block[f"{k}/fit_seconds"].mean()), 3),
                              cindex_multimodal=round(float(block[f"{k}/multimodal"].mean()), 2))
                         for m, k in (("CCA-MLE", "cca-mle"), ("EM", "em"))])
    cost.to_csv(os.path.join(args.out, "final_nll_and_cost.csv"), index=False)

    pd.set_option("display.width", 200)
    print(kpi[kpi.setting == "multimodal"].to_string(index=False))
    print()
    print(tt.to_string(index=False))
    print()
    print(cost.to_string(index=False))
    print(f"\nseeds averaged per fold: {block['n_seeds']}; wrote {args.out}")


if __name__ == "__main__":
    main()
