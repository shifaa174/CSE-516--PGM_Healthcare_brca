"""
make_figures.py -- Result figures for the M2 report.
Run AFTER run_experiments.py and sweep_latent_dim.py.

    python make_figures.py ../../Results
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch

S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, SURF = "#0b0b0b", "#52514e", "#b8b7b2", "#ffffff"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF,
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.labelcolor": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.bbox": "tight",
})


def fig1_graphical_model(out):
    fig, ax = plt.subplots(figsize=(5.6, 3.1))
    ax.set_xlim(0.2, 10); ax.set_ylim(0, 6); ax.axis("off")
    ax.add_patch(Rectangle((0.7, 0.5), 6.5, 4.6, fill=False, ec=INK2, lw=1.0))
    ax.text(6.95, 0.78, r"$i = 1 \ldots N$", ha="right", va="bottom",
            fontsize=8, color=INK2)

    def node(x, y, label, shaded, color):
        ax.add_patch(Circle((x, y), 0.62,
                            facecolor=color if shaded else "white",
                            edgecolor=INK, lw=1.3, zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=11,
                color="white" if shaded else INK, zorder=4)

    node(2.1, 3.6, r"$z_i$", False, "white")
    node(2.1, 1.6, r"$x_i$", True, S1)
    node(4.6, 1.6, r"$y_i$", True, S2)
    node(6.0, 3.6, r"$h_i$", False, "white")
    node(8.8, 3.6, r"$T_i,\delta_i$", True, S3)

    for a, b in [((2.1, 3.6), (2.1, 1.6)), ((2.1, 3.6), (4.6, 1.6)),
                 ((2.1, 3.6), (6.0, 3.6)), ((6.0, 3.6), (8.8, 3.6))]:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12,
                                     shrinkA=20, shrinkB=20, lw=1.2,
                                     color=INK, zorder=2))

    ax.text(1.78, 2.6, r"$W_x,\Psi_x$", fontsize=8, color=INK2, ha="right")
    ax.text(3.9, 2.75, r"$W_y,\Psi_y$", fontsize=8, color=INK2)
    ax.text(4.05, 3.85, r"$\beta$-free", fontsize=8, color=INK2, ha="center")
    ax.text(7.4, 3.85, r"Cox PH $\beta$", fontsize=8, color=INK2, ha="center")
    ax.text(2.1, 4.55, "latent", ha="center", fontsize=7.5, color=INK2)
    ax.text(2.1, 0.72, "genomics", ha="center", fontsize=7.5, color=INK2)
    ax.text(4.6, 0.72, "histology", ha="center", fontsize=7.5, color=INK2)
    ax.text(8.8, 2.65, "survival", ha="center", fontsize=7.5, color=INK2)

    ax.set_title("Implemented graphical model (Group G02)\n"
                 r"shaded $=$ observed, unshaded $=$ latent/derived;  "
                 r"$X \perp\!\!\!\perp Y \mid Z$",
                 fontsize=9, color=INK, loc="left")
    fig.savefig(os.path.join(out, "fig1_graphical_model.png"))
    plt.close(fig)


def fig2_em_convergence(res, out):
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    traces = [r["em"]["loglik_trace"] for r in res["per_fold"]]
    L = max(len(t) for t in traces)
    for t in traces:
        ax.plot(np.arange(1, len(t) + 1), np.array(t) / 1e3,
                color=MUTED, lw=1.0, zorder=2)
    padded = np.empty((len(traces), L))
    for i, t in enumerate(traces):
        padded[i, :len(t)] = t
        padded[i, len(t):] = t[-1]
    ax.plot(np.arange(1, L + 1), padded.mean(0) / 1e3, color=S1, lw=2.0,
            zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel("EM iteration")
    ax.set_ylabel(r"training log-likelihood ($\times 10^{3}$ nats)")
    ax.grid(axis="y", color=MUTED, lw=0.5, alpha=0.5); ax.set_axisbelow(True)
    mono = all(r["em"]["loglik_monotone"] for r in res["per_fold"])
    iters = [r["em"]["n_iter"] for r in res["per_fold"]]
    ax.set_title("EM increases the joint log-likelihood monotonically\n"
                 f"bold = mean of {len(traces)} folds; grey = per fold; "
                 f"monotone on every fold: {mono}; "
                 f"{min(iters)}–{max(iters)} iterations",
                 fontsize=9, loc="left", color=INK)
    fig.savefig(os.path.join(out, "fig2_em_convergence.png"))
    plt.close(fig)


def fig3_cindex(res, out):
    s = res["summary"]
    settings = ["multimodal", "genomics_only", "histology_only"]
    labels = ["Multimodal\n$p(z|x,y)$", "Genomics only\n$p(z|x)$",
              "Histology only\n$p(z|y)$"]
    base = ["baseline_pca_concat", "baseline_pca_genomics",
            "baseline_pca_histology"]
    series = [("CCA-MLE (baseline route)", S1,
               [s[f"c_index/cca-mle/{k}"] for k in settings]),
              ("EM (our extension)", S2,
               [s[f"c_index/em/{k}"] for k in settings]),
              ("PCA + Cox (non-PGM baseline)", S3,
               [s[f"c_index/{k}"] for k in base])]

    fig, ax = plt.subplots(figsize=(6.3, 3.5))
    x = np.arange(len(settings)); w = 0.26
    for i, (name, col, vals) in enumerate(series):
        mu = [v["mean"] for v in vals]
        sd = [v["std"] for v in vals]
        pos = x + (i - 1) * (w + 0.015)
        ax.bar(pos, mu, w, label=name, color=col, edgecolor=SURF, lw=1.2,
               zorder=3)
        ax.errorbar(pos, mu, yerr=sd, fmt="none", ecolor=INK2, elinewidth=1.0,
                    capsize=3, zorder=4)
        for p_, m_ in zip(pos, mu):
            ax.text(p_, m_ + 1.3, f"{m_:.1f}", ha="center", va="bottom",
                    fontsize=7.5, color=INK2, zorder=5)
    ax.axhline(50, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(-0.44, 50.9, "random (50)", fontsize=7.5, color=INK2, ha="left")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("test C-index (%)"); ax.set_ylim(40, 86)
    ax.grid(axis="y", color=MUTED, lw=0.5, alpha=0.5); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.16))
    cfg = res["config"]
    ax.set_title("Survival prediction from the latent posterior\n"
                 f"mean ± sd over {len(cfg['folds'])} folds, "
                 f"variant = {cfg['variant']}, d = {cfg['d']}, "
                 f"seed = {cfg['seed']}",
                 fontsize=9, loc="left", color=INK, pad=30)
    fig.savefig(os.path.join(out, "fig3_cindex_comparison.png"))
    plt.close(fig)


def fig4_latent_dim(sweep, out):
    """Two panels, never a dual axis: the two measures have different units."""
    grid = sorted({r["d"] for r in sweep["rows"]})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for meth, col, lab in [("cca-mle", S1, "CCA-MLE"), ("em", S2, "EM")]:
        ll = [sweep["summary"][f"{meth}/d={d}"]["valid_loglik_mean"]
              for d in grid]
        ci = [sweep["summary"][f"{meth}/d={d}"]["c_index_mean"] for d in grid]
        sd = [sweep["summary"][f"{meth}/d={d}"]["c_index_std"] for d in grid]
        a1.plot(grid, ll, "o-", color=col, lw=2, ms=5, label=lab, zorder=3)
        a2.errorbar(grid, ci, yerr=sd, fmt="o-", color=col, lw=2, ms=5,
                    ecolor=INK2, elinewidth=0.9, capsize=3, label=lab, zorder=3)
    for ax, yl in ((a1, "validation log-likelihood / sample"),
                   (a2, "test C-index (%)")):
        ax.set_xlabel("latent dimension $d$"); ax.set_ylabel(yl)
        ax.grid(color=MUTED, lw=0.5, alpha=0.5); ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8)
    a1.set_title("Unsupervised criterion: monotone in $d$", fontsize=9,
                 loc="left", color=INK)
    a2.set_title("Supervised outcome: peaks near $d=10$", fontsize=9,
                 loc="left", color=INK)
    fig.suptitle("Held-out likelihood does not identify a useful latent "
                 "dimension", fontsize=10, x=0.02, ha="left", y=1.06,
                 color=INK)
    fig.savefig(os.path.join(out, "fig4_latent_dim_sweep.png"))
    plt.close(fig)


def fig5_entropy(res, out):
    s = res["summary"]
    settings = ["multimodal", "genomics_only", "histology_only"]
    labels = ["Multimodal", "Genomics only", "Histology only"]
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    x = np.arange(len(settings)); w = 0.36
    for i, (meth, col, lab) in enumerate([("cca-mle", S1, "CCA-MLE"),
                                          ("em", S2, "EM")]):
        v = [s[f"posterior_entropy/{meth}/{k}"]["mean"] for k in settings]
        pos = x + (i - 0.5) * (w + 0.02)
        ax.bar(pos, v, w, color=col, edgecolor=SURF, lw=1.2, label=lab,
               zorder=3)
        for p_, m_ in zip(pos, v):
            ax.text(p_, m_ + 0.15, f"{m_:.2f}", ha="center", va="bottom",
                    fontsize=7.5, color=INK2, zorder=4)
    d = res["config"]["d"]
    prior = 0.5 * d * np.log(2 * np.pi * np.e)
    ax.axhline(prior, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(len(settings) - 0.55, prior + 0.2,
            f"prior entropy ({prior:.1f})", fontsize=7.5, color=INK2,
            ha="right")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("posterior differential entropy (nats)")
    ax.grid(axis="y", color=MUTED, lw=0.5, alpha=0.5); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.set_title("Latent uncertainty rises when a modality is withheld\n"
                 r"$H[p(z\mid\cdot)]$, mean over folds",
                 fontsize=9, loc="left", color=INK)
    fig.savefig(os.path.join(out, "fig5_posterior_entropy.png"))
    plt.close(fig)


if __name__ == "__main__":
    out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "../../Results")
    res = json.load(open(os.path.join(out, "m2_results_main.json")))
    fig1_graphical_model(out)
    fig2_em_convergence(res, out)
    fig3_cindex(res, out)
    fig5_entropy(res, out)
    sp = os.path.join(out, "m2_latent_dim_sweep.json")
    if os.path.exists(sp):
        fig4_latent_dim(json.load(open(sp)), out)
    print("figures written to", out)
