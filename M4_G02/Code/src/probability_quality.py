import argparse, os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from sksurv.metrics import brier_score, integrated_brier_score
from sksurv.util import Surv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data import load_fold, standardise
from pcca import ProbabilisticCCA

warnings.filterwarnings("ignore")


def features(method, tr, te, d, seed, noise):
    if method == "Non-PGM PCA baseline":
        return np.hstack([tr["X"], tr["Y"]]), np.hstack([te["X"], te["Y"]])
    m = ProbabilisticCCA(d=d, reg=1e-2, random_state=seed, noise=noise)
    if method.startswith("CCA-MLE"):
        m.fit_cca_mle(tr["X"], tr["Y"])
    else:
        m.fit_em(tr["X"], tr["Y"], max_iter=300, tol=1e-7, init="random")
    return m.map_estimate(X=tr["X"], Y=tr["Y"]), m.map_estimate(X=te["X"], Y=te["Y"])


def frame(Z, times, events):
    df = pd.DataFrame(Z, columns=[f"z{i}" for i in range(Z.shape[1])])
    df["T"] = np.maximum(np.asarray(times, float), 1.0)
    df["E"] = np.asarray(events, int)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--variant", default="paperpca")
    ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--d-cca", type=int, default=30)
    ap.add_argument("--d-em", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--penalizer", type=float, default=0.1)
    ap.add_argument("--horizons", type=float, nargs="+", default=[1095, 1826])
    ap.add_argument("--bins", type=int, default=5)
    ap.add_argument("--out", default="M4_G02/Results/probability_quality")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    methods = [("CCA-MLE", args.d_cca, "block"), ("EM", args.d_em, "block"),
               ("CCA-MLE diagonal", args.d_cca, "diag"), ("Non-PGM PCA baseline", None, None)]
    rows, calib = [], []
    for fold in args.folds:
        raw = load_fold(args.root, fold, args.variant)
        tr, _, te = standardise(raw["train"], raw["valid"], raw["test"])
        for name, d, noise in methods:
            Ztr, Zte = features(name, tr, te, d, args.seed + fold, noise)
            ftr, fte = frame(Ztr, tr["times"], tr["events"]), frame(Zte, te["times"], te["events"])
            cph = CoxPHFitter(penalizer=args.penalizer).fit(ftr, "T", "E")
            ytr = Surv.from_arrays(ftr["E"].astype(bool), ftr["T"])
            yte = Surv.from_arrays(fte["E"].astype(bool), fte["T"])
            hz = [h for h in args.horizons if h < fte["T"].max() and h < ftr["T"].max()]
            sf = cph.predict_survival_function(fte, times=hz).T.values
            _, bs = brier_score(ytr, yte, sf, hz)
            grid = np.linspace(np.percentile(fte["T"], 10), np.percentile(fte["T"], 80), 30)
            ibs = integrated_brier_score(ytr, yte, cph.predict_survival_function(fte, times=grid).T.values, grid)
            row = dict(fold=fold, method=name,
                       c_index=100 * cph.score(fte, scoring_method="concordance_index"),
                       test_partial_nll_per_sample=-cph.score(fte, scoring_method="log_likelihood"),
                       integrated_brier=ibs)
            for h, b in zip(hz, bs):
                row[f"brier_{int(h)}d"] = b
            rows.append(row)
            for j, h in enumerate(hz):
                p = sf[:, j]
                q = pd.qcut(p, args.bins, labels=False, duplicates="drop")
                for b in np.unique(q):
                    idx = q == b
                    km = KaplanMeierFitter().fit(fte["T"][idx], fte["E"][idx])
                    calib.append(dict(fold=fold, method=name, horizon_days=int(h), bin=int(b),
                                      n=int(idx.sum()), predicted=float(p[idx].mean()),
                                      observed=float(km.predict(h))))
        print(f"fold {fold} done", flush=True)

    per_fold = pd.DataFrame(rows)
    per_fold.to_csv(os.path.join(args.out, "probability_quality_per_fold.csv"), index=False)
    cal = pd.DataFrame(calib)
    cal["abs_gap"] = (cal["predicted"] - cal["observed"]).abs()
    cal.to_csv(os.path.join(args.out, "calibration_bins.csv"), index=False)
    ece = (cal.assign(w=cal["abs_gap"] * cal["n"]).groupby(["method", "horizon_days"])
           .apply(lambda g: g["w"].sum() / g["n"].sum()).rename("ece").reset_index())
    ece.to_csv(os.path.join(args.out, "calibration_ece.csv"), index=False)

    metrics = [c for c in per_fold.columns if c not in ("fold", "method")]
    agg = per_fold.groupby("method", sort=False)[metrics].agg(["mean", "std"]).round(4)
    agg.columns = [f"{a}_{b}" for a, b in agg.columns]
    agg.to_csv(os.path.join(args.out, "probability_quality_summary.csv"))

    hs = sorted(cal["horizon_days"].unique())
    fig, axes = plt.subplots(1, len(hs), figsize=(5.5 * len(hs), 4.8))
    axes = np.atleast_1d(axes)
    for ax, h in zip(axes, hs):
        ax.plot([0, 1], [0, 1], "--", color="grey", lw=1, label="perfect calibration")
        for name, _, _ in methods:
            c = cal[(cal.method == name) & (cal.horizon_days == h)].groupby("bin")[["predicted", "observed"]].mean()
            ax.plot(c["predicted"], c["observed"], "o-", label=name)
        ax.set_title(f"Calibration at {round(h / 365.25)} years")
        ax.set_xlabel("predicted survival probability")
        ax.set_ylabel("observed survival (Kaplan-Meier)")
        ax.set_xlim(0.4, 1.02)
        ax.set_ylim(0.4, 1.02)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "calibration_plot.png"), dpi=150)

    pd.set_option("display.width", 200)
    print(agg.to_string())
    print(ece.to_string(index=False))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
