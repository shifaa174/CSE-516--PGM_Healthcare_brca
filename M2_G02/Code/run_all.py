import os, subprocess, sys

if len(sys.argv) < 2:
    sys.exit("usage: python run_all.py <path to repo containing 03_brca/>")

root = os.path.abspath(sys.argv[1])
here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, "src")
out = os.path.join(os.path.dirname(here), "Results")
py = sys.executable

steps = [
    ("1/6  paper baseline: authors' trained models and PCA input",
     ["paper_baseline.py", "--root", root, "--out", out]),
    ("2/6  latent-dimension sweep, paper PCA input",
     ["sweep_latent_dim.py", "--root", root, "--variant", "paperpca", "--grid", "5", "10", "20", "30", "50",
      "--em-iters", "120", "--seed", "42", "--out", out, "--tag", "paperpca"]),
    ("3/6  our PGM on the paper PCA input, d chosen on validation",
     ["run_experiments.py", "--root", root, "--variant", "paperpca", "--noise", "block",
      "--d-from", os.path.join(out, "m2_latent_dim_sweep_paperpca.json"), "--seed", "42", "--out", out, "--tag", "paperpca"]),
    ("4/6  latent-dimension sweep, reduced input",
     ["sweep_latent_dim.py", "--root", root, "--variant", "pca", "--grid", "2", "5", "10", "20", "30",
      "--em-iters", "120", "--seed", "42", "--out", out, "--tag", "main"]),
    ("5/6  our PGM on the reduced input, d chosen on validation",
     ["run_experiments.py", "--root", root, "--variant", "pca", "--noise", "block",
      "--d-from", os.path.join(out, "m2_latent_dim_sweep_main.json"), "--seed", "42", "--out", out, "--tag", "main"]),
    ("6/6  figures and summary", ["make_figures.py", out]),
    (None, ["write_summary.py", out]),
]

os.makedirs(out, exist_ok=True)
for title, args in steps:
    if title:
        print(f"== {title} ==", flush=True)
    subprocess.run([py] + args, cwd=src, check=True)
print(f"done -- results in {out}")
