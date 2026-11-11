import argparse, json, os, platform, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data import load_fold, standardise
from pcca import ProbabilisticCCA
from survival import CoxPH, c_index

SETTINGS = {"multimodal":     ("X", "Y"),
            "genomics_only":  ("X",),
            "histology_only": ("Y",)}


def downstream_cindex(Ztr, tr, Zte, te, penalizer):
    cox = CoxPH(penalizer=penalizer).fit(Ztr, tr["times"], tr["events"])
    return float(c_index(te["times"], te["events"], cox.risk(Zte)))


def select_latent_dim(tr, va, grid, args):

    scores = {}
    for d in grid:
        m = ProbabilisticCCA(d=d, reg=args.reg, random_state=args.seed,
                             noise=args.noise).fit_cca_mle(tr["X"], tr["Y"])
        scores[d] = m.score_per_sample(va["X"], va["Y"])
    best = max(scores, key=scores.get)
    return best, {int(k): float(v) for k, v in scores.items()}


def run_fold(args, fold):
    raw = load_fold(args.root, fold, args.variant)
    tr, va, te = standardise(raw["train"], raw["valid"], raw["test"])

    res = {"fold": fold,
           "n_train": int(len(tr["times"])),
           "n_valid": int(len(va["times"])),
           "n_test":  int(len(te["times"])),
           "p_genomic":   int(tr["X"].shape[1]),
           "p_histology": int(tr["Y"].shape[1]),
           "events_train": int(tr["events"].sum()),
           "events_test":  int(te["events"].sum()),
           "event_rate_train": float(tr["events"].mean()),
           "event_rate_test":  float(te["events"].mean())}

    dims = {"cca-mle": args.d, "em": args.d}
    if args.d_from:
        dims.update({k: int(v) for k, v in json.load(open(args.d_from, encoding="utf-8"))["chosen_d"].items()})
    if args.select_d:
        d, curve = select_latent_dim(tr, va, args.d_grid, args)
        dims = {"cca-mle": d, "em": d}
        res["latent_dim_selection"] = {"chosen": int(d),
                                       "val_loglik_per_sample": curve}
    res["latent_dim"] = {k: int(v) for k, v in dims.items()}

    for method in ("cca-mle", "em"):
        d = dims[method]
        m = ProbabilisticCCA(d=d, reg=args.reg, random_state=args.seed + fold,
                             noise=args.noise)
        t0 = time.time()
        if method == "cca-mle":
            m.fit_cca_mle(tr["X"], tr["Y"])
        else:
            m.fit_em(tr["X"], tr["Y"], max_iter=args.em_iters, tol=args.em_tol,
                     init=args.em_init, verbose=args.verbose)
        fit_s = time.time() - t0

        trace = [float(x) for x in m.loglik_trace_]
        monotone = all(b >= a - 1e-6 for a, b in zip(trace, trace[1:]))

        blk = {"latent_dim": int(d),
               "fit_seconds": round(fit_s, 3),
               "n_iter": int(m.n_iter_),
               "loglik_monotone": bool(monotone),
               "train_loglik_per_sample": m.score_per_sample(tr["X"], tr["Y"]),
               "valid_loglik_per_sample": m.score_per_sample(va["X"], va["Y"]),
               "test_loglik_per_sample":  m.score_per_sample(te["X"], te["Y"]),
               "loglik_trace": trace}
        if method == "cca-mle":
            blk["canonical_correlations"] = [float(r) for r in m.rho_]

        t0 = time.time()
        for name, keys in SETTINGS.items():
            kw_tr = {k: tr[k] for k in keys}
            kw_te = {k: te[k] for k in keys}
            Ztr, _ = m.posterior(**kw_tr)
            Zte, G = m.posterior(**kw_te)
            _, logdet = np.linalg.slogdet(G)
            blk[name] = {
                "c_index": downstream_cindex(Ztr, tr, Zte, te, args.penalizer),
                "posterior_entropy_nats":
                    float(0.5 * (d * np.log(2 * np.pi * np.e) + logdet)),
            }
        blk["inference_seconds"] = round(time.time() - t0, 3)
        res[method] = blk

    res["baseline_pca_concat"] = downstream_cindex(
        np.hstack([tr["X"], tr["Y"]]), tr,
        np.hstack([te["X"], te["Y"]]), te, args.penalizer)
    res["baseline_pca_genomics"] = downstream_cindex(
        tr["X"], tr, te["X"], te, args.penalizer)
    res["baseline_pca_histology"] = downstream_cindex(
        tr["Y"], tr, te["Y"], te, args.penalizer)
    return res


def summarise(all_res):
    def agg(getter, scale=100.0):
        v = [getter(r) for r in all_res]
        return {"mean": round(float(np.mean(v)) * scale, 3),
                "std":  round(float(np.std(v)) * scale, 3),
                "per_fold": [round(float(x) * scale, 3) for x in v]}

    s = {}
    for meth in ("cca-mle", "em"):
        for setting in SETTINGS:
            s[f"c_index/{meth}/{setting}"] = agg(
                lambda r, m=meth, t=setting: r[m][t]["c_index"])
        s[f"loglik_per_sample/{meth}/test"] = agg(
            lambda r, m=meth: r[m]["test_loglik_per_sample"], scale=1.0)
        s[f"loglik_per_sample/{meth}/train"] = agg(
            lambda r, m=meth: r[m]["train_loglik_per_sample"], scale=1.0)
        s[f"fit_seconds/{meth}"] = agg(
            lambda r, m=meth: r[m]["fit_seconds"], scale=1.0)
        s[f"n_iter/{meth}"] = agg(lambda r, m=meth: r[m]["n_iter"], scale=1.0)
        for setting in SETTINGS:
            s[f"posterior_entropy/{meth}/{setting}"] = agg(
                lambda r, m=meth, t=setting:
                    r[m][t]["posterior_entropy_nats"], scale=1.0)
    for b in ("baseline_pca_concat", "baseline_pca_genomics",
              "baseline_pca_histology"):
        s[f"c_index/{b}"] = agg(lambda r, k=b: r[k])
    return s


def main():
    ap = argparse.ArgumentParser(description="M2 pipeline, Group G02")
    ap.add_argument("--root", required=True,
                    help="repository root containing 03_brca/")
    ap.add_argument("--variant", default="pca", choices=["pca", "histogram", "paperpca", "loghist"])
    ap.add_argument("--noise", default="block", choices=["block", "diag"])
    ap.add_argument("--out", default="../../Results")
    ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--d", type=int, default=10, help="latent dimension")
    ap.add_argument("--d-from", default=None,
                    help="sweep JSON whose validation-chosen d is used per route")
    ap.add_argument("--select-d", action="store_true",
                    help="select d on the validation split")
    ap.add_argument("--d-grid", type=int, nargs="+",
                    default=[2, 5, 10, 20, 30])
    ap.add_argument("--reg", type=float, default=1e-2)
    ap.add_argument("--penalizer", type=float, default=0.1)
    ap.add_argument("--em-iters", type=int, default=300)
    ap.add_argument("--em-tol", type=float, default=1e-7)
    ap.add_argument("--em-init", default="random", choices=["random", "cca"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", default=None, help="suffix for the output file")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    np.random.seed(args.seed)
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)

    all_res = []
    for f in args.folds:
        print(f"[fold {f}] running ...", flush=True)
        all_res.append(run_fold(args, f))

    payload = {
        "config": {k: v for k, v in vars(args).items()},
        "environment": {"python": platform.python_version(),
                        "numpy": np.__version__,
                        "platform": platform.platform()},
        "per_fold": all_res,
        "summary": summarise(all_res),
    }
    tag = args.tag or f"{args.variant}_{args.noise}"
    path = os.path.join(out, f"m2_results_{tag}.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)

    print("\n=== C-index (mean +/- sd over folds, %) ===")
    for k, v in payload["summary"].items():
        if k.startswith("c_index/"):
            print(f"  {k[8:]:38s} {v['mean']:6.2f} +/- {v['std']:.2f}")
    print("\n=== Held-out log-likelihood per sample ===")
    for m in ("cca-mle", "em"):
        v = payload["summary"][f"loglik_per_sample/{m}/test"]
        print(f"  {m:8s} {v['mean']:10.2f} +/- {v['std']:.2f}")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
