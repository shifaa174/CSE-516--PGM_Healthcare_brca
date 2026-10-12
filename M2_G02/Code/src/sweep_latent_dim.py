"""
sweep_latent_dim.py -- Representation KPI: effect of the latent dimension d.
===============================================================================
d is the project's only structural hyperparameter. The M1 plan was to select it
on the validation split by held-out log-likelihood. This script tests whether
that criterion actually identifies a useful d, by sweeping d and recording BOTH
the unsupervised criterion (validation log-likelihood) and the supervised
outcome (test C-index) for each learning route.

Usage
    python sweep_latent_dim.py --root <repo> --out ../../Results
"""
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data import load_fold, standardise        # noqa: E402
from pcca import ProbabilisticCCA              # noqa: E402
from survival import CoxPH, c_index            # noqa: E402


def cindex_from(m, tr, te, penalizer):
    Ztr, _ = m.posterior(X=tr["X"], Y=tr["Y"])
    Zte, _ = m.posterior(X=te["X"], Y=te["Y"])
    cox = CoxPH(penalizer=penalizer).fit(Ztr, tr["times"], tr["events"])
    return float(c_index(te["times"], te["events"], cox.risk(Zte)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--variant", default="pca", choices=["pca", "histogram"])
    ap.add_argument("--noise", default="block", choices=["block", "diag"])
    ap.add_argument("--grid", type=int, nargs="+", default=[2, 5, 10, 20, 30, 50])
    ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--penalizer", type=float, default=0.1)
    ap.add_argument("--em-iters", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="../../Results")
    args = ap.parse_args()

    rows = []
    for fold in args.folds:
        raw = load_fold(args.root, fold, args.variant)
        tr, va, te = standardise(raw["train"], raw["valid"], raw["test"])
        for d in args.grid:
            for method in ("cca-mle", "em"):
                m = ProbabilisticCCA(d=d, random_state=args.seed + fold,
                                     noise=args.noise)
                if method == "cca-mle":
                    m.fit_cca_mle(tr["X"], tr["Y"])
                else:
                    m.fit_em(tr["X"], tr["Y"], max_iter=args.em_iters)
                rows.append({
                    "fold": fold, "d": d, "method": method,
                    "valid_loglik_per_sample": m.score_per_sample(va["X"], va["Y"]),
                    "test_loglik_per_sample":  m.score_per_sample(te["X"], te["Y"]),
                    "test_c_index": cindex_from(m, tr, te, args.penalizer),
                })
            print(f"  fold {fold}, d={d} done", flush=True)

    summary = {}
    for method in ("cca-mle", "em"):
        for d in args.grid:
            sel = [r for r in rows if r["method"] == method and r["d"] == d]
            summary[f"{method}/d={d}"] = {
                "valid_loglik_mean": round(float(np.mean(
                    [r["valid_loglik_per_sample"] for r in sel])), 3),
                "c_index_mean": round(float(np.mean(
                    [r["test_c_index"] for r in sel])) * 100, 2),
                "c_index_std": round(float(np.std(
                    [r["test_c_index"] for r in sel])) * 100, 2),
            }

    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "m2_latent_dim_sweep.json")
    with open(path, "w") as fh:
        json.dump({"config": vars(args), "rows": rows, "summary": summary},
                  fh, indent=2)

    print("\n=== Latent dimension sweep (mean over folds) ===")
    print(f"{'route':9s} {'d':>4} {'valid logL/n':>13} {'test C-index':>14}")
    for k, v in summary.items():
        meth, dd = k.split("/d=")
        print(f"{meth:9s} {dd:>4} {v['valid_loglik_mean']:>13.2f} "
              f"{v['c_index_mean']:>8.2f} +/- {v['c_index_std']:.2f}")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
