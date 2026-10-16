import argparse, contextlib, io, json, os, platform, sys, warnings
import numpy as np
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from survival import CoxPH, c_index

ARXIV_TABLE4 = {
    "PCA": {"genomics": (65.76, 9.02), "imaging": (61.24, 6.39), "concat": (66.86, 10.40)},
    "SCCA/HD": {"genomics": (65.95, 10.70), "imaging": (54.00, 11.23), "concat": (65.71, 7.96)},
    "SCCA/PD": {"genomics": (68.34, 8.40), "imaging": (55.66, 11.14), "concat": (67.37, 6.58)},
    "SCCA/OPD": {"genomics": (67.63, 9.43), "imaging": (54.19, 11.23), "concat": (65.25, 5.61)},
    "GN-SCCA/PD": {"genomics": (64.94, 8.98), "imaging": (56.81, 13.31), "concat": (66.42, 12.08)},
    "GN-SCCA/OPD": {"genomics": (65.52, 7.67), "imaging": (57.87, 10.28), "concat": (68.08, 9.88)},
}

CONFIGS = [("PCA", None, None)] + [
    (f"{lab}/{d.upper()}", m, d)
    for lab, m in (("SCCA", "scca"), ("GN-SCCA", "gnscca"))
    for d in ("hd", "pd", "opd")
]
SETTINGS = ("genomics", "imaging", "concat")


def stored_npz(brca, label, method, defl, setting):
    sr = os.path.join(brca, "survival_results")
    if label == "PCA":
        name = {"genomics": "PCAgenomics", "imaging": "PCAimaging", "concat": "PCAconcat"}[setting]
    else:
        name = f"feature_{method}_log_histogram_{defl}" + ("" if setting == "concat" else f"_{setting}")
    p = os.path.join(sr, name + ".npz")
    if not os.path.exists(p):
        return None
    try:
        a = np.load(p, allow_pickle=True)["arr_0"].astype(float) * 100
        return a if np.all(np.isfinite(a)) else None
    except Exception:
        return None


def paper_engine_available():
    try:
        import lifelines, sksurv, pandas, sklearn
        return True
    except Exception:
        return False


def run_with_paper_code(brca):
    code_dir = os.path.join(brca, "A_survival_prediction")
    here = os.getcwd()
    sys.path.insert(0, code_dir)
    os.chdir(code_dir)
    try:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import run_coxph_final as R
        mode_map = {"PCA": {"genomics": "PCAgenomics", "imaging": "PCAimaging", "concat": "PCAconcat"}}
        out = {}
        for label, method, defl in CONFIGS:
            for setting in SETTINGS:
                if label == "PCA":
                    st = {"mode": mode_map["PCA"][setting], "method": "scca",
                          "suffix": "log_histogram", "deflation": "opd"}
                else:
                    st = {"mode": "feature" if setting == "concat" else setting,
                          "method": method, "suffix": "log_histogram", "deflation": defl}
                try:
                    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        scores, _ = R.run_all_fold(st, range(5))
                    arr = np.array([np.nan if x is None else x for x in scores], float) * 100
                except Exception:
                    arr = np.full(5, np.nan)
                out[(label, setting)] = arr
        return out
    finally:
        os.chdir(here)
        sys.path.remove(code_dir)


def _load(brca, split, fold):
    base = os.path.join(brca, "data", "splits_journal")
    lh = os.path.join(base, f"{split}_log_histogram", "1000", str(fold))
    pc = os.path.join(base, f"{split}_pca", "100", str(fold))
    rd = lambda p: np.loadtxt(p, delimiter=",", ndmin=2)
    return dict(G=rd(os.path.join(lh, "genes.csv")), I=rd(os.path.join(lh, "image.csv")),
                PG=rd(os.path.join(pc, "genes.csv")), PI=rd(os.path.join(lh, "PCA_image.csv")),
                t=np.loadtxt(os.path.join(lh, "times.csv"), ndmin=1),
                e=np.loadtxt(os.path.join(lh, "events.csv"), ndmin=1))


def run_with_our_code(brca, penalizer):
    out = {}
    folds = [(_load(brca, "train", k), _load(brca, "test", k)) for k in range(5)]
    for label, method, defl in CONFIGS:
        for setting in SETTINGS:
            scores = []
            for k, (tr, te) in enumerate(folds):
                try:
                    if label == "PCA":
                        Ftr = {"genomics": tr["PG"], "imaging": tr["PI"], "concat": np.hstack([tr["PG"], tr["PI"]])}[setting]
                        Fte = {"genomics": te["PG"], "imaging": te["PI"], "concat": np.hstack([te["PG"], te["PI"]])}[setting]
                    else:
                        prefix = "sparse" if method == "scca" else method
                        m = loadmat(os.path.join(brca, "results", f"{prefix}_{defl}_1000_{k}_log_histogram_best.mat"))
                        U, V = np.squeeze(m["u"]), np.squeeze(m["v"])
                        xa, ya = tr["G"] @ U, tr["I"] @ V
                        xb, yb = te["G"] @ U, te["I"] @ V
                        Ftr = {"genomics": xa, "imaging": ya, "concat": np.hstack([xa, ya])}[setting]
                        Fte = {"genomics": xb, "imaging": yb, "concat": np.hstack([xb, yb])}[setting]
                    if not (np.all(np.isfinite(Ftr)) and np.all(np.isfinite(Fte))):
                        raise ValueError
                    cox = CoxPH(penalizer=penalizer).fit(Ftr, tr["t"], tr["e"])
                    scores.append(c_index(te["t"], te["e"], cox.risk(Fte)) * 100)
                except Exception:
                    scores.append(np.nan)
            out[(label, setting)] = np.array(scores, float)
    return out


def ms(a):
    if a is None or not np.all(np.isfinite(a)):
        return None
    return {"mean": round(float(a.mean()), 2), "std": round(float(a.std()), 2),
            "per_fold": [round(float(x), 2) for x in a]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--engine", default="auto", choices=["auto", "paper", "ours"])
    ap.add_argument("--penalizer", type=float, default=0.1)
    ap.add_argument("--out", default="../../Results")
    args = ap.parse_args()

    brca = os.path.join(os.path.abspath(args.root), "03_brca")
    engine = args.engine
    if engine == "auto":
        engine = "paper" if paper_engine_available() else "ours"
    if engine == "paper" and not paper_engine_available():
        sys.exit("paper engine needs: pip3 install lifelines scikit-survival pandas scikit-learn")
    print(f"engine: {engine} ({'authors survival code, lifelines Cox' if engine == 'paper' else 'our Cox on the authors trained models'})")

    res = run_with_paper_code(brca) if engine == "paper" else run_with_our_code(brca, args.penalizer)

    rows = []
    for label, method, defl in CONFIGS:
        for setting in SETTINGS:
            ours = ms(res[(label, setting)])
            stored = ms(stored_npz(brca, label, method, defl, setting))
            ax = ARXIV_TABLE4.get(label, {}).get(setting)
            rows.append({"method": label, "setting": setting, "rerun": ours, "authors_stored": stored,
                         "arxiv_table4": {"mean": ax[0], "std": ax[1]} if ax else None})

    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    payload = {"engine": engine, "config": vars(args),
               "environment": {"python": platform.python_version(), "numpy": np.__version__,
                               "platform": platform.platform()},
               "rows": rows}
    path = os.path.join(out, "m2_paper_baseline.json")
    json.dump(payload, open(path, "w", encoding="utf-8"), indent=2)

    f = lambda d: "    failed    " if d is None else f"{d['mean']:6.2f} ± {d['std']:5.2f}"
    print(f"\n{'method':12s} {'setting':9s} {'rerun':>15s} {'authors stored':>16s} {'arXiv Table 4':>15s}")
    for r in rows:
        print(f"{r['method']:12s} {r['setting']:9s} {f(r['rerun']):>15s} {f(r['authors_stored']):>16s} {f(r['arxiv_table4']):>15s}")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
