#!/usr/bin/env python3
"""Run additional missing-modality, observation-noise and training-size experiments.

Run from repository root:
  python M3_Analysis_Code/run_robustness.py --root . --variant paperpca --out M3_G02/Results/robustness

Uses the existing M2 data.py, pcca.py and survival.py without modifying them.
"""
import argparse, csv, json, sys, time
from pathlib import Path
import numpy as np

def load_modules(root):
    src = Path(root) / "M2_G02" / "Code" / "src"
    if not src.exists():
        raise FileNotFoundError(f"Could not find M2_G02/Code/src at {src.resolve()}")
    sys.path.insert(0, str(src.resolve()))
    from data import load_fold, standardise
    from pcca import ProbabilisticCCA
    from survival import CoxPH, c_index
    return load_fold, standardise, ProbabilisticCCA, CoxPH, c_index

def posterior_rows(model, X, Y, mode):
    """Infer one patient at a time so each row can have a different observed modality."""
    z = []
    for i in range(len(X)):
        if mode == "both":
            zi = model.map_estimate(X=X[i:i+1], Y=Y[i:i+1])
        elif mode == "genomics_only":
            zi = model.map_estimate(X=X[i:i+1])
        elif mode == "histology_only":
            zi = model.map_estimate(Y=Y[i:i+1])
        else:
            raise ValueError(mode)
        z.append(zi[0])
    return np.asarray(z)

def fit_model(method, X, Y, d, seed, noise, reg, em_iters):
    from pcca import ProbabilisticCCA
    model = ProbabilisticCCA(d=d, reg=reg, random_state=seed, noise=noise)
    if method == "cca-mle":
        model.fit_cca_mle(X, Y)
    else:
        model.fit_em(X, Y, max_iter=em_iters, tol=1e-7, init="random")
    return model

def score(model, cox, X, Y, times, events, mode, c_index):
    Z = posterior_rows(model, X, Y, mode)
    return float(c_index(times, events, cox.risk(Z)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--variant", default="paperpca",
                    choices=["paperpca","pca","histogram","loghist"])
    ap.add_argument("--folds", nargs="+", type=int, default=[0,1,2,3,4])
    ap.add_argument("--methods", nargs="+", default=["cca-mle","em"],
                    choices=["cca-mle","em"])
    ap.add_argument("--d-cca", type=int, default=30)
    ap.add_argument("--d-em", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--noise-model", choices=["block","diag"], default="block")
    ap.add_argument("--reg", type=float, default=1e-2)
    ap.add_argument("--cox-penalizer", type=float, default=0.1)
    ap.add_argument("--em-iters", type=int, default=300)
    ap.add_argument("--missing-rates", nargs="+", type=float, default=[0.0,0.1,0.25,0.5])
    ap.add_argument("--noise-levels", nargs="+", type=float, default=[0.0,0.1,0.25,0.5,1.0])
    ap.add_argument("--train-fractions", nargs="+", type=float, default=[0.25,0.5,0.75,1.0])
    ap.add_argument("--out", default="M3_G02/Results/robustness")
    args = ap.parse_args()

    load_fold, standardise, ProbabilisticCCA, CoxPH, c_index = load_modules(args.root)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    rng_master = np.random.default_rng(args.seed)

    for fold in args.folds:
        raw = load_fold(args.root, fold, args.variant)
        tr, va, te = standardise(raw["train"], raw["valid"], raw["test"])
        for method in args.methods:
            d = args.d_cca if method == "cca-mle" else args.d_em
            model = ProbabilisticCCA(d=d, reg=args.reg,
                                     random_state=args.seed+fold,
                                     noise=args.noise_model)
            start = time.time()
            if method == "cca-mle":
                model.fit_cca_mle(tr["X"], tr["Y"])
            else:
                model.fit_em(tr["X"], tr["Y"], max_iter=args.em_iters,
                             tol=1e-7, init="random")
            fit_seconds = time.time()-start
            Ztr = model.map_estimate(X=tr["X"], Y=tr["Y"])
            cox = CoxPH(penalizer=args.cox_penalizer).fit(
                Ztr, tr["times"], tr["events"])

            # Reference: complete paired modalities.
            ref = score(model, cox, te["X"], te["Y"], te["times"],
                        te["events"], "both", c_index)
            rows.append(dict(fold=fold, method=method, experiment="reference",
                             level=0.0, c_index=ref, fit_seconds=fit_seconds,
                             n_train=len(tr["times"]), d=d))

            # Missing-modality robustness. For each rate, randomly choose patients
            # whose genomic OR histology modality is hidden, balanced by random draw.
            for rate in args.missing_rates:
                rng = np.random.default_rng(args.seed + 1000*fold + int(rate*10000))
                n = len(te["times"])
                nmask = int(round(rate*n))
                mask = np.zeros(n, dtype=bool)
                if nmask:
                    mask[rng.choice(n, size=nmask, replace=False)] = True
                z = []
                for i in range(n):
                    if not mask[i]:
                        zi = model.map_estimate(X=te["X"][i:i+1], Y=te["Y"][i:i+1])
                    elif rng.random() < 0.5:
                        zi = model.map_estimate(X=te["X"][i:i+1])
                    else:
                        zi = model.map_estimate(Y=te["Y"][i:i+1])
                    z.append(zi[0])
                risk = cox.risk(np.asarray(z))
                ci = float(c_index(te["times"], te["events"], risk))
                rows.append(dict(fold=fold, method=method,
                                 experiment="random_one_modality_missing",
                                 level=rate, c_index=ci, fit_seconds=fit_seconds,
                                 n_train=len(tr["times"]), d=d))

            # Observation-noise robustness. Standardised inputs have unit-ish SD.
            for level in args.noise_levels:
                rng = np.random.default_rng(args.seed + 2000*fold + int(level*10000))
                Xn = te["X"] + rng.normal(0, level, te["X"].shape)
                Yn = te["Y"] + rng.normal(0, level, te["Y"].shape)
                ci = score(model, cox, Xn, Yn, te["times"], te["events"],
                           "both", c_index)
                rows.append(dict(fold=fold, method=method,
                                 experiment="gaussian_observation_noise",
                                 level=level, c_index=ci, fit_seconds=fit_seconds,
                                 n_train=len(tr["times"]), d=d))

            # Training-size sensitivity. Keep the original held-out test fold fixed.
            for fraction in args.train_fractions:
                if not (0 < fraction <= 1):
                    raise ValueError("Training fractions must be in (0,1].")
                rng = np.random.default_rng(args.seed + 3000*fold + int(fraction*10000))
                nsub = max(10, int(round(fraction*len(tr["times"]))))
                ix = rng.choice(len(tr["times"]), size=min(nsub,len(tr["times"])),
                                replace=False)
                sub = {k: tr[k][ix] for k in ("X","Y","times","events")}
                # Guard against latent dimensions exceeding available sample rank.
                dsub = min(d, sub["X"].shape[1], sub["Y"].shape[1], len(ix)-1)
                msub = ProbabilisticCCA(d=dsub, reg=args.reg,
                                        random_state=args.seed+fold,
                                        noise=args.noise_model)
                if method == "cca-mle":
                    msub.fit_cca_mle(sub["X"], sub["Y"])
                else:
                    msub.fit_em(sub["X"], sub["Y"], max_iter=args.em_iters,
                                tol=1e-7, init="random")
                zsub = msub.map_estimate(X=sub["X"], Y=sub["Y"])
                coxsub = CoxPH(penalizer=args.cox_penalizer).fit(
                    zsub, sub["times"], sub["events"])
                ci = score(msub, coxsub, te["X"], te["Y"], te["times"],
                           te["events"], "both", c_index)
                rows.append(dict(fold=fold, method=method,
                                 experiment="training_fraction",
                                 level=fraction, c_index=ci,
                                 fit_seconds=np.nan, n_train=len(ix), d=dsub))
            print(f"Fold {fold}, {method}: completed", flush=True)

    # Write raw rows and mean +/- SD across fixed folds for each condition.
    raw_path = out / "robustness_per_fold.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    groups = {}
    for r in rows:
        key = (r["method"], r["experiment"], r["level"])
        groups.setdefault(key, []).append(r["c_index"])
    summary = []
    for (method, experiment, level), vals in sorted(groups.items()):
        a = np.asarray(vals, float)
        summary.append(dict(method=method, experiment=experiment, level=level,
                            n_folds=len(a), c_index_mean=float(np.nanmean(a)),
                            c_index_sd=float(np.nanstd(a, ddof=1)) if len(a)>1 else 0.0,
                            c_index_mean_percent=float(100*np.nanmean(a)),
                            c_index_sd_percentage_points=float(100*np.nanstd(a, ddof=1)) if len(a)>1 else 0.0))
    with (out/"robustness_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader(); w.writerows(summary)
    (out/"robustness_config.json").write_text(json.dumps(vars(args), indent=2), encoding="utf-8")
    print(f"Saved raw and summary robustness outputs to {out.resolve()}")

if __name__ == "__main__":
    main()
