import os, subprocess, sys

if len(sys.argv) < 2:
    sys.exit("usage: python run_all.py <path to repo containing 03_brca/>")

root = os.path.abspath(sys.argv[1])
here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, "src")
out = os.path.join(os.path.dirname(here), "Results")
py = sys.executable
seeds = ["42", "43", "44", "45", "46"]
sweep = os.path.join(out, "latent_dim_sweep_paperpca.json")
ablations = [("abl_noise_diag", ["--noise", "diag"]), ("abl_eminit_cca", ["--em-init", "cca"]),
             ("abl_reg_1e-3", ["--reg", "1e-3"]), ("abl_reg_1e-1", ["--reg", "1e-1"]),
             ("abl_pen_0.01", ["--penalizer", "0.01"]), ("abl_pen_1.0", ["--penalizer", "1.0"])]

steps = [
    ("1/10  paper baselines (PCA, SCCA, GN-SCCA)",
     [["paper_baseline.py", "--root", root, "--out", out]]),
    ("2/10  latent-dimension sweep, paper PCA input",
     [["sweep_latent_dim.py", "--root", root, "--variant", "paperpca", "--grid", "5", "10", "20", "30", "50",
       "--em-iters", "120", "--seed", "42", "--out", out, "--tag", "paperpca"]]),
    ("3/10  latent-dimension sweep and PGM, reduced input",
     [["sweep_latent_dim.py", "--root", root, "--variant", "pca", "--grid", "2", "5", "10", "20", "30",
       "--em-iters", "120", "--seed", "42", "--out", out, "--tag", "main"],
      ["run_experiments.py", "--root", root, "--variant", "pca", "--noise", "block",
       "--d-from", os.path.join(out, "latent_dim_sweep_main.json"), "--seed", "42", "--out", out, "--tag", "main"]]),
    ("4/10  main results, paper PCA input, 5 seeds x 5 folds",
     [["run_experiments.py", "--root", root, "--variant", "paperpca", "--noise", "block", "--d-from", sweep,
       "--seed", s, "--out", out, "--tag", f"paperpca_seed{s}"] for s in seeds]
     + [["run_experiments.py", "--root", root, "--variant", "paperpca", "--noise", "block", "--d-from", sweep,
         "--seed", "42", "--out", out, "--tag", "paperpca"]]),
    ("5/10  ablations, 5 seeds each",
     [["run_experiments.py", "--root", root, "--variant", "paperpca", "--d-from", sweep,
       "--seed", s, "--out", out, "--tag", f"{tag}_seed{s}"] + flag for tag, flag in ablations for s in seeds]),
    ("6/10  robustness: missing modality, noise, training size",
     [["run_robustness.py", "--root", root, "--out", os.path.join(out, "robustness")]]),
    ("7/10  aggregation and exploratory tests",
     [["aggregate_results.py", "--results", out, "--out", os.path.join(out, "analysis")]]),
    ("8/10  final statistics: 95% CI, paired t-test, Wilcoxon",
     [["final_stats.py", "--results", out, "--out", os.path.join(out, "stats")]]),
    ("9/10  probability quality: NLL, Brier, calibration",
     [["probability_quality.py", "--root", root, "--out", os.path.join(out, "probability_quality")]]),
    ("10/10 figures and summary",
     [["make_figures.py", out], ["write_summary.py", out],
      ["plot_results.py", "--analysis", os.path.join(out, "analysis"),
       "--robustness", os.path.join(out, "robustness"), "--out", os.path.join(out, "plots")]]),
]

os.makedirs(out, exist_ok=True)
for title, cmds in steps:
    print(f"== {title} ==", flush=True)
    for args in cmds:
        subprocess.run([py] + args, cwd=src, check=True, stdout=subprocess.DEVNULL)
print(f"done -- results in {out}")
