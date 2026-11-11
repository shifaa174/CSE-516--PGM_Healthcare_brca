import os
import numpy as np

VARIANTS = {
    "pca": (("pca", "100", "genes.csv"), ("pca", "100", "image.csv")),
    "histogram": (("histogram", "1000", "genes.csv"), ("histogram", "1000", "image.csv")),
    "paperpca": (("pca", "100", "genes.csv"), ("log_histogram", "1000", "PCA_image.csv")),
    "loghist": (("log_histogram", "1000", "genes.csv"), ("log_histogram", "1000", "image.csv")),
}


def _dir(root, split, suffix, size, fold):
    return os.path.join(root, "03_brca", "data", "splits_journal",
                        f"{split}_{suffix}", size, str(fold))


def load_split(root, split, fold=0, variant="histogram"):
    (xs, xz, xf), (ys, yz, yf) = VARIANTS[variant]
    px = _dir(root, split, xs, xz, fold)
    py = _dir(root, split, ys, yz, fold)
    for p in (px, py):
        if not os.path.isdir(p):
            raise FileNotFoundError(
                f"split directory not found: {p}\n"
                f"Check --root points at the repository containing 03_brca/.")
    return dict(
        X=np.loadtxt(os.path.join(px, xf), delimiter=",", ndmin=2),
        Y=np.loadtxt(os.path.join(py, yf), delimiter=",", ndmin=2),
        times=np.loadtxt(os.path.join(py, "times.csv"), ndmin=1),
        events=np.loadtxt(os.path.join(py, "events.csv"), ndmin=1),
    )


def load_fold(root, fold=0, variant="histogram"):
    return {s: load_split(root, s, fold, variant)
            for s in ("train", "valid", "test")}


def standardise(train, *others):

    stats = {k: (train[k].mean(0), train[k].std(0) + 1e-8) for k in ("X", "Y")}

    def apply(d):
        out = dict(d)
        for k in ("X", "Y"):
            mu, sd = stats[k]
            out[k] = (d[k] - mu) / sd
        return out

    return [apply(train)] + [apply(o) for o in others]
