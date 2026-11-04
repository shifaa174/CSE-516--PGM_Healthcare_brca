#!/usr/bin/env python3
"""Create report-ready plots from aggregate and robustness CSV outputs."""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--analysis", default="M3_G02/Results/analysis")
    ap.add_argument("--robustness", default="M3_G02/Results/robustness")
    ap.add_argument("--out", default="M3_G02/Results/plots")
    a=ap.parse_args()
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    analysis=Path(a.analysis)
    cpath=analysis/"cindex_summary_percent.csv"
    if cpath.exists():
        df=pd.read_csv(cpath)
        df=df[df.metric.str.contains("multimodal|genomics_only|histology_only",regex=True)]
        for outcome in ("multimodal","genomics_only","histology_only"):
            sub=df[df.metric.str.endswith(outcome)]
            if sub.empty: continue
            fig,ax=plt.subplots(figsize=(10,5))
            x=range(len(sub))
            means=sub.mean_pct.to_numpy()
            errs=sub.mean_fold_sd.to_numpy()*100
            ax.bar(list(x),means,yerr=errs,capsize=4)
            ax.set_xticks(list(x))
            ax.set_xticklabels(sub.experiment.tolist(),rotation=35,ha="right")
            ax.set_ylabel("C-index (%)")
            ax.set_title(f"M3 C-index comparison: {outcome.replace('_',' ')}")
            ax.grid(axis="y",alpha=.25)
            fig.tight_layout()
            fig.savefig(out/f"cindex_{outcome}.png",dpi=220)
            plt.close(fig)

    rp=Path(a.robustness)/"robustness_summary.csv"
    if rp.exists():
        r=pd.read_csv(rp)
        for experiment,title,xlabel in [
            ("random_one_modality_missing","Missing-modality robustness","Fraction of test patients with one modality hidden"),
            ("gaussian_observation_noise","Observation-noise robustness","Gaussian noise SD added to standardised features"),
            ("training_fraction","Training-size sensitivity","Fraction of training patients used"),
        ]:
            sub=r[r.experiment==experiment]
            if sub.empty: continue
            fig,ax=plt.subplots(figsize=(7,4.5))
            for method,g in sub.groupby("method"):
                g=g.sort_values("level")
                ax.errorbar(g.level,g.c_index_mean_percent,
                            yerr=g.c_index_sd_percentage_points,
                            marker="o",capsize=3,label=method)
            ax.set_xlabel(xlabel); ax.set_ylabel("Test C-index (%)")
            ax.set_title(title); ax.legend(); ax.grid(alpha=.25)
            fig.tight_layout(); fig.savefig(out/f"{experiment}.png",dpi=220)
            plt.close(fig)
    print(f"Plots saved in {out.resolve()}")

if __name__=="__main__":
    main()
