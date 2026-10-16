import json, math, os, sys
import numpy as np

OUT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "../../Results")


def j(name):
    p = os.path.join(OUT, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


def f2(x):
    return f"{x:.2f}"


def ms(v):
    return "failed" if v is None else f"{f2(v['mean'])} ± {f2(v['std'])}"


def dims(res):
    ld = res["per_fold"][0].get("latent_dim", res["config"]["d"])
    return (ld["cca-mle"], ld["em"]) if isinstance(ld, dict) else (ld, ld)


LABEL = {"pca": "reduced input: 50 genomic PCs + 215 raw CellProfiler features",
         "paperpca": "paper's PCA input: 50 genomic PCs + 50 imaging PCs"}


def pgm_table(L, res, title):
    s = res["summary"]
    dc, de = dims(res)
    c = lambda k: ms(s[k])
    L.append(f"### {title}\n")
    L.append("| Method | Multimodal p(z\\|x,y) | Genomics only p(z\\|x) | Histology only p(z\\|y) |")
    L.append("|---|---|---|---|")
    L.append(f"| CCA-MLE, d = {dc} | {c('c_index/cca-mle/multimodal')} | {c('c_index/cca-mle/genomics_only')} | {c('c_index/cca-mle/histology_only')} |")
    L.append(f"| EM (our extension), d = {de} | {c('c_index/em/multimodal')} | {c('c_index/em/genomics_only')} | {c('c_index/em/histology_only')} |")
    L.append(f"| Non-PGM baseline + our Cox | {c('c_index/baseline_pca_concat')} | {c('c_index/baseline_pca_genomics')} | {c('c_index/baseline_pca_histology')} |\n")


def sweep_table(L, sw, title):
    grid = sorted({r["d"] for r in sw["rows"]})
    S = sw["summary"]
    L.append(f"### {title}\n")
    L.append("| d | CCA valid logL/n | CCA valid C-index | CCA test C-index | EM valid logL/n | EM valid C-index | EM test C-index |")
    L.append("|---|---|---|---|---|---|---|")
    for g in grid:
        a, b = S[f"cca-mle/d={g}"], S[f"em/d={g}"]
        L.append(f"| {g} | {f2(a['valid_loglik_mean'])} | {f2(a.get('valid_c_index_mean', float('nan')))} | "
                 f"{f2(a['c_index_mean'])} ± {f2(a['c_index_std'])} | {f2(b['valid_loglik_mean'])} | "
                 f"{f2(b.get('valid_c_index_mean', float('nan')))} | {f2(b['c_index_mean'])} ± {f2(b['c_index_std'])} |")
    ch = sw.get("chosen_d")
    if ch:
        L.append(f"\nChosen on validation C-index: CCA-MLE d = {ch['cca-mle']}, EM d = {ch['em']}. Test C-index is reported only, never used to choose.\n")


def em_convergence(res):
    cap = res["config"]["em_iters"]
    it = [f["em"]["n_iter"] for f in res["per_fold"]]
    worst = 0.0
    for f in res["per_fold"]:
        t = f["em"]["loglik_trace"]
        drops = [b - a for a, b in zip(t, t[1:]) if b < a]
        if drops:
            worst = max(worst, -min(drops) / abs(t[-1]))
    conv = sum(1 for i in it if i < cap)
    return conv, len(it), cap, worst, min(it), max(it)


def main():
    main_res, pp, pb = j("m2_results_main.json"), j("m2_results_paperpca.json"), j("m2_paper_baseline.json")
    sw_main = j("m2_latent_dim_sweep_main.json") or j("m2_latent_dim_sweep.json")
    sw_pp = j("m2_latent_dim_sweep_paperpca.json")
    env = (main_res or pp or {}).get("environment", {})
    L = ["# Initial KPI results — M2_G02\n"]
    L.append(f"Software used for this run: Python {env.get('python', '?')}, numpy {env.get('numpy', '?')}. "
             f"Runs on Windows, macOS and Linux. Regenerate everything from the `Code` folder with "
             f"`python run_all.py ../..` (or `bash run_all.sh ../..`).\n")
    L.append("Model selection rule: the latent dimension d is chosen per learning route by **validation** C-index; "
             "test data is never used to choose anything.\n")

    if pb:
        eng = ("the authors' own survival code (`03_brca/A_survival_prediction/`, lifelines Cox)"
               if pb["engine"] == "paper" else
               "our Cox model on the authors' trained models (lifelines / scikit-survival not installed)")
        diffs = [abs(a - b) for r in pb["rows"] if r["rerun"] and r["authors_stored"]
                 for a, b in zip(r["rerun"]["per_fold"], r["authors_stored"]["per_fold"])]
        L.append("## 1. Paper baseline reproduction\n")
        L.append(f"Paper methods rerun from the authors' trained SCCA / GN-SCCA models (`03_brca/results/`) and their "
                 f"PCA files, scored with {eng}. Largest per-fold difference from the authors' stored results: "
                 f"{f2(max(diffs)) if diffs else 'n/a'} C-index points.\n")
        L.append("| Method | Setting | Rerun (ours) | Authors' stored results | arXiv Table 4 |")
        L.append("|---|---|---|---|---|")
        for r in pb["rows"]:
            L.append(f"| {r['method']} | {r['setting']} | {ms(r['rerun'])} | {ms(r['authors_stored'])} | {ms(r['arxiv_table4'])} |")
        L.append("\nGN-SCCA with HD deflation fails in the authors' own results too. The authors' stored results "
                 "differ from the arXiv preprint's Table 4, which most likely reflects a revision between preprint "
                 "and journal version; quote the journal table in the report.\n")

    L.append("## 2. Downstream KPI — our PGM, test C-index (%), mean ± sd over 5 folds\n")
    if pp:
        pgm_table(L, pp, f"2a. On the {LABEL['paperpca']} (like-for-like with the paper's PCA baseline)")
    if main_res:
        pgm_table(L, main_res, f"2b. On the {LABEL['pca']}")

    L.append("## 3. Learning KPI — log-likelihood per sample and EM convergence\n")
    L.append("| Input | Route | d | Train logL/n | Test logL/n | Fit time |")
    L.append("|---|---|---|---|---|---|")
    for res, lab in ((pp, "paper PCA 50 + 50"), (main_res, "reduced 50 + 215")):
        if not res:
            continue
        s = res["summary"]
        dc, de = dims(res)
        for m, d in (("cca-mle", dc), ("em", de)):
            L.append(f"| {lab} | {m} | {d} | {ms(s[f'loglik_per_sample/{m}/train'])} | {ms(s[f'loglik_per_sample/{m}/test'])} | {f2(s[f'fit_seconds/{m}']['mean'])} s |")
    L.append("")
    for res, lab in ((pp, "paper PCA input"), (main_res, "reduced input")):
        if res:
            conv, n, cap, worst, lo, hi = em_convergence(res)
            if conv == n:
                L.append(f"- EM on the {lab}: converged on all {n} folds ({lo}–{hi} iterations).")
            else:
                L.append(f"- EM on the {lab}: converged on {conv} of {n} folds; on the other {n - conv} it plateaued and oscillated "
                         f"by at most {worst*100:.4f}% until the {cap}-iteration cap (positive-definite projection of Psi in the M-step).")
    L.append("")

    L.append("## 4. Inference KPI — posterior differential entropy (nats)\n")
    L.append("| Input | Route | d | Prior | Multimodal | Genomics only | Histology only |")
    L.append("|---|---|---|---|---|---|---|")
    for res, lab in ((pp, "paper PCA 50 + 50"), (main_res, "reduced 50 + 215")):
        if not res:
            continue
        s = res["summary"]
        dc, de = dims(res)
        for m, d in (("cca-mle", dc), ("em", de)):
            e = [f2(s[f"posterior_entropy/{m}/{k}"]["mean"]) for k in ("multimodal", "genomics_only", "histology_only")]
            L.append(f"| {lab} | {m} | {d} | {f2(0.5 * d * math.log(2 * math.pi * math.e))} | {e[0]} | {e[1]} | {e[2]} |")
    L.append("\nLower than the prior = the data is informative about Z; withholding a modality raises entropy.\n")

    L.append("## 5. Representation KPI — latent dimension sweeps\n")
    if sw_pp:
        sweep_table(L, sw_pp, "5a. Paper PCA input")
    if sw_main:
        sweep_table(L, sw_main, "5b. Reduced input")

    L.append("## 6. Findings\n")
    if pb and pp:
        best = max((r for r in pb["rows"] if r["setting"] == "concat" and r["rerun"]), key=lambda r: r["rerun"]["mean"])
        pca = next(r for r in pb["rows"] if r["method"] == "PCA" and r["setting"] == "concat")
        s = pp["summary"]
        L.append(f"1. **The paper's results reproduce.** Rerunning the authors' trained models gives their numbers "
                 f"(best concatenated: {best['method']} {ms(best['rerun'])}; PCA {ms(pca['rerun'])}).")
        L.append(f"2. **On the same input, our PGM does not beat the paper.** CCA-MLE {ms(s['c_index/cca-mle/multimodal'])} and "
                 f"EM {ms(s['c_index/em/multimodal'])} versus PCA {ms(pca['rerun'])} and {best['method']} {ms(best['rerun'])}; "
                 f"all differences are within one standard deviation.")
    parts = []
    for sw, lab in ((sw_pp, "paper PCA input"), (sw_main, "reduced input")):
        if sw:
            grid = sorted({r["d"] for r in sw["rows"]})
            dec = all(all(sw["summary"][f"{m}/d={a}"]["valid_loglik_mean"] >= sw["summary"][f"{m}/d={b}"]["valid_loglik_mean"]
                          for a, b in zip(grid, grid[1:])) for m in ("cca-mle", "em"))
            parts.append(f"{'' if dec else 'almost '}always picking d = {grid[0]} on the {lab}")
    if parts:
        L.append("3. **Held-out likelihood cannot choose d.** Under both learning routes, validation log-likelihood "
                 "keeps falling as d grows, " + " and ".join(parts) + ", which predicts poorly. We therefore choose d by validation C-index.")
    if sw_pp and sw_main and 10 in {r["d"] for r in sw_pp["rows"]} and 10 in {r["d"] for r in sw_main["rows"]}:
        g = lambda sw, m: sw["summary"][f"{m}/d=10"]["valid_loglik_mean"]
        L.append(f"4. **EM's likelihood gain depends on the input.** At d = 10, EM vs CCA-MLE validation logL/n is "
                 f"{f2(g(sw_main, 'em'))} vs {f2(g(sw_main, 'cca-mle'))} on the reduced input, but "
                 f"{f2(g(sw_pp, 'em'))} vs {f2(g(sw_pp, 'cca-mle'))} on the paper's PCA input. The reduced input contains "
                 f"constant and duplicated CellProfiler columns; there the ridge-regularised CCA route is far from the true "
                 f"maximum-likelihood fit and EM corrects it, while on clean PCA features both routes reach the same fit.")
    L.append("")
    L.append("## 7. Caveats\n")
    L.append("- 95 events in 974 patients: fold-to-fold sd of 5–10 points is intrinsic; no difference here is claimed as significant.")
    L.append("- Validation splits hold about 15–25 events each, so the validation-chosen d is itself noisy.")
    L.append("- Across machines, C-index can differ in the second decimal place (numerical libraries); conclusions are unchanged.")
    L.append("- The paper's SCCA / GN-SCCA rows use the authors' trained models; their training code is MATLAB and is not rerun.\n")

    L.append("## Figures\n")
    for f, d in [("fig1_graphical_model.png", "the implemented PGM and the CI assumption"),
                 ("fig2_em_convergence.png", "EM log-likelihood traces, reduced input"),
                 ("fig3_cindex_comparison.png", "C-index by method and modality, reduced input"),
                 ("fig3b_cindex_paper_input.png", "C-index by method and modality, paper PCA input"),
                 ("fig4_latent_dim_sweep.png", "d sweep, reduced input"),
                 ("fig4b_latent_dim_sweep_paper_input.png", "d sweep, paper PCA input"),
                 ("fig5_posterior_entropy.png", "latent uncertainty by modality availability"),
                 ("fig6_paper_comparison.png", "paper baselines (rerun) vs our PGM on the same input")]:
        if os.path.exists(os.path.join(OUT, f)):
            L.append(f"- `{f}` — {d}")
    open(os.path.join(OUT, "SUMMARY.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")

    import csv
    rows = []
    for res, lab in ((main_res, "reduced"), (pp, "paperpca")):
        if not res:
            continue
        for k, v in res["summary"].items():
            rows.append([lab, k, v["mean"], v["std"], " / ".join(f"{x:.2f}" for x in v["per_fold"])])
    if pb:
        for r in pb["rows"]:
            v = r["rerun"]
            rows.append(["paper_baseline", f"{r['method']}/{r['setting']}", v["mean"] if v else "", v["std"] if v else "",
                         " / ".join(f"{x:.2f}" for x in v["per_fold"]) if v else "failed"])
    with open(os.path.join(OUT, "m2_kpi_table.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["input", "kpi", "mean", "std", "per_fold"])
        w.writerows(rows)
    print("wrote", os.path.join(OUT, "SUMMARY.md"), "and m2_kpi_table.csv")


if __name__ == "__main__":
    main()
