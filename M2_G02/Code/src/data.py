"""
data.py -- Dataset pipeline.
===============================================================================
Loads the processed TCGA-BRCA splits released with the base paper.

Expected layout, relative to --root:

    03_brca/data/splits_journal/
        {train,valid,test}_histogram/1000/{0..4}/   full features
        {train,valid,test}_pca/100/{0..4}/          reduced features
            genes.csv    genomic feature matrix   X
            image.csv    histology feature matrix Y
            times.csv    survival time in days    T
            events.csv   event indicator (1=death, 0=censored)  delta
            patients.csv TCGA patient barcodes

Variants
    "histogram"  1000 genes x 1075 CellProfiler histogram features  (primary)
    "pca"        50 genomic PCs x 215 histology PCs                 (fast)

LEAKAGE CONTROL: `standardise` fits the z-scoring statistics on the TRAINING
split only and applies them unchanged to validation and test. Patients are
never shared across splits -- the splits ship that way from the base paper.
"""
import os
import numpy as np

VARIANTS = {"pca": ("pca", "100"), "histogram": ("histogram", "1000")}


def _dir(root, split, variant, fold):
    suffix, size = VARIANTS[variant]
    return os.path.join(root, "03_brca", "data", "splits_journal",
                        f"{split}_{suffix}", size, str(fold))


def load_split(root, split, fold=0, variant="histogram"):
    p = _dir(root, split, variant, fold)
    if not os.path.isdir(p):
        raise FileNotFoundError(
            f"split directory not found: {p}\n"
            f"Check --root points at the repository containing 03_brca/.")
    return dict(
        X=np.loadtxt(os.path.join(p, "genes.csv"), delimiter=",", ndmin=2),
        Y=np.loadtxt(os.path.join(p, "image.csv"), delimiter=",", ndmin=2),
        times=np.loadtxt(os.path.join(p, "times.csv"), ndmin=1),
        events=np.loadtxt(os.path.join(p, "events.csv"), ndmin=1),
    )


def load_fold(root, fold=0, variant="histogram"):
    return {s: load_split(root, s, fold, variant)
            for s in ("train", "valid", "test")}


def standardise(train, *others):
    """Fit z-scoring on TRAIN ONLY; apply to every split. No leakage."""
    stats = {k: (train[k].mean(0), train[k].std(0) + 1e-8) for k in ("X", "Y")}

    def apply(d):
        out = dict(d)
        for k in ("X", "Y"):
            mu, sd = stats[k]
            out[k] = (d[k] - mu) / sd
        return out

    return [apply(train)] + [apply(o) for o in others]
