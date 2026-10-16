import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data import load_fold, standardise
from pcca import ProbabilisticCCA
from survival import CoxPH, c_index


def cindex_from(m, tr, te, penalizer, va=None):
    Ztr, _ = m.posterior(X=tr["X"], Y=tr["Y"])
    Zte, _ = m.posterior(X=te["X"], Y=te["Y"])
    cox = CoxPH(penalizer=penalizer).fit(Ztr, tr["times"], tr["events"])
    test = float(c_index(te["times"], te["events"], cox.risk(Zte)))
    if va is None:
        return test
    Zva, _ = m.posterior(X=va["X"], Y=va["Y"])
    return test, float(c_index(va["times"], va["events"], cox.risk(Zva)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--variant", default="pca", choices=["pca", "histogram", "paperpca", "loghist"])
    ap.add_argument("--noise", default="block", choices=["block", "diag"])
    ap.add_argument("--grid", type=int, nargs="+", default=[2, 5, 10, 20, 30, 50])
    ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--penalizer", type=float, default=0.1)
    ap.add_argument("--em-iters", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="../../Results")
    ap.add_argument("--tag", default=None)
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
                    "valid_c_index": cindex_from(m, tr, te, args.penalizer, va)[1],
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
                "valid_c_index_mean": round(float(np.mean(
                    [r["valid_c_index"] for r in sel])) * 100, 2),
            }

    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    chosen = {}
    for meth in ("cca-mle", "em"):
        cand = {d: summary[f"{meth}/d={d}"]["valid_c_index_mean"] for d in args.grid}
        chosen[meth] = int(max(cand, key=cand.get))
    name = "m2_latent_dim_sweep" + (f"_{args.tag}" if args.tag else "") + ".json"
    path = os.path.join(out, name)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"config": vars(args), "rows": rows, "summary": summary,
                   "chosen_d": chosen, "selection": "validation C-index"},
                  fh, indent=2)

    print("\n=== Latent dimension sweep (mean over folds) ===")
    print(f"{'route':9s} {'d':>4} {'valid logL/n':>13} {'valid C-idx':>12} {'test C-index':>14}")
    for k, v in summary.items():
        meth, dd = k.split("/d=")
        print(f"{meth:9s} {dd:>4} {v['valid_loglik_mean']:>13.2f} {v['valid_c_index_mean']:>12.2f} "
              f"{v['c_index_mean']:>8.2f} +/- {v['c_index_std']:.2f}")

    for meth in ("cca-mle", "em"):
        print(f"  {meth}: d chosen by validation C-index = {chosen[meth]}")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
