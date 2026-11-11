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
    iters = [r["em"]["n_iter"] for r in res["per_fold"]]
    cap = res["config"]["em_iters"]
    conv = sum(1 for i in iters if i < cap)
    ax.set_title("EM log-likelihood rises steeply, then plateaus\n"
                 f"bold = mean of {len(traces)} folds; grey = per fold; "
                 f"converged on {conv} of {len(iters)} folds before the {cap}-iteration cap",
                 fontsize=9, loc="left", color=INK)
    fig.savefig(os.path.join(out, "fig2_em_convergence.png"))
    plt.close(fig)


def _dims(res):
    ld = res["per_fold"][0].get("latent_dim", res["config"]["d"])
    if isinstance(ld, dict):
        return ld["cca-mle"], ld["em"]
    return ld, ld


def _variant_label(v):
    return {"pca": "50 genomic PCs + 215 raw CellProfiler features",
            "paperpca": "paper's 50 + 50 PCA input",
            "histogram": "1000 genes + 1075 histology features",
            "loghist": "1000 genes + 1075 log-histogram features"}.get(v, v)


def fig3_cindex(res, out, name="fig3_cindex_comparison.png"):
    s = res["summary"]
    settings = ["multimodal", "genomics_only", "histology_only"]
    labels = ["Multimodal\n$p(z|x,y)$", "Genomics only\n$p(z|x)$",
              "Histology only\n$p(z|y)$"]
    base = ["baseline_pca_concat", "baseline_pca_genomics",
            "baseline_pca_histology"]
    dc, de = _dims(res)
    series = [(f"CCA-MLE, d = {dc}", S1,
               [s[f"c_index/cca-mle/{k}"] for k in settings]),
              (f"EM (our extension), d = {de}", S2,
               [s[f"c_index/em/{k}"] for k in settings]),
              ("PCA + Cox (non-PGM baseline)", S3,
               [s[f"c_index/{k}"] for k in base])]

    fig, ax = plt.subplots(figsize=(6.3, 3.5))
    x = np.arange(len(settings)); w = 0.26
    for i, (slabel, col, vals) in enumerate(series):
        mu = [v["mean"] for v in vals]
        sd = [v["std"] for v in vals]
        pos = x + (i - 1) * (w + 0.015)
        ax.bar(pos, mu, w, label=slabel, color=col, edgecolor=SURF, lw=1.2,
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
                 f"mean ± sd over {len(cfg['folds'])} folds; input: {_variant_label(cfg['variant'])}; "
                 f"seed {cfg['seed']}",
                 fontsize=9, loc="left", color=INK, pad=30)
    fig.savefig(os.path.join(out, name))
    plt.close(fig)


def fig4_latent_dim(sweep, out, name="fig4_latent_dim_sweep.png"):
    grid = sorted({r["d"] for r in sweep["rows"]})
    chosen = sweep.get("chosen_d", {})
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.0))
    panels = [("valid_loglik_mean", None, "validation log-likelihood / sample", "Likelihood: falls as d grows"),
              ("valid_c_index_mean", None, "validation C-index (%)", "Validation C-index: picks d"),
              ("c_index_mean", "c_index_std", "test C-index (%)", "Test C-index: reported only")]
    for ax, (key, sdk, yl, title) in zip(axes, panels):
        for meth, col, lab in [("cca-mle", S1, "CCA-MLE"), ("em", S2, "EM")]:
            ys = [sweep["summary"][f"{meth}/d={d}"].get(key, np.nan) for d in grid]
            if sdk:
                sd = [sweep["summary"][f"{meth}/d={d}"][sdk] for d in grid]
                ax.errorbar(grid, ys, yerr=sd, fmt="o-", color=col, lw=2, ms=5,
                            ecolor=INK2, elinewidth=0.9, capsize=3, label=lab, zorder=3)
            else:
                ax.plot(grid, ys, "o-", color=col, lw=2, ms=5, label=lab, zorder=3)
            if key == "valid_c_index_mean" and meth in chosen:
                d0 = chosen[meth]
                y0 = sweep["summary"][f"{meth}/d={d0}"][key]
                ax.scatter([d0], [y0], s=140, facecolor="none", edgecolor=col, lw=2, zorder=4)
        ax.set_xlabel("latent dimension $d$"); ax.set_ylabel(yl)
        ax.grid(color=MUTED, lw=0.5, alpha=0.5); ax.set_axisbelow(True)
        ax.set_title(title, fontsize=9, loc="left", color=INK)
    axes[0].legend(frameon=False, fontsize=8)
    v = sweep["config"].get("variant", "pca")
    ch = ", ".join(f"{k} d = {v2}" for k, v2 in chosen.items())
    fig.suptitle(f"Latent dimension sweep, input: {_variant_label(v)}. Circled = chosen on validation ({ch})",
                 fontsize=10, x=0.01, ha="left", y=1.06, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, name))
    plt.close(fig)


def fig5_entropy(res, out):
    s = res["summary"]
    settings = ["multimodal", "genomics_only", "histology_only"]
    labels = ["Multimodal", "Genomics only", "Histology only"]
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    x = np.arange(len(settings)); w = 0.36
    for i, (meth, col, lab) in enumerate([("cca-mle", S1, f"CCA-MLE, d = {_dims(res)[0]}"),
                                          ("em", S2, f"EM, d = {_dims(res)[1]}")]):
        v = [s[f"posterior_entropy/{meth}/{k}"]["mean"] for k in settings]
        pos = x + (i - 0.5) * (w + 0.02)
        ax.bar(pos, v, w, color=col, edgecolor=SURF, lw=1.2, label=lab,
               zorder=3)
        for p_, m_ in zip(pos, v):
            ax.text(p_, m_ + 0.15, f"{m_:.2f}", ha="center", va="bottom",
                    fontsize=7.5, color=INK2, zorder=4)
    dc, de = _dims(res)
    if dc == de:
        prior = 0.5 * dc * np.log(2 * np.pi * np.e)
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


def fig6_paper_comparison(pb, pp, out):
    pick = {r["method"]: r for r in pb["rows"] if r["setting"] == "concat"}
    items = []
    for lab in ("PCA", "SCCA/OPD", "GN-SCCA/OPD"):
        r = pick.get(lab)
        if r and r["rerun"]:
            items.append((f"Paper: {lab}", S3, r["rerun"]["mean"], r["rerun"]["std"]))
    if pp:
        s = pp["summary"]
        dc, de = _dims(pp)
        items.append((f"Ours: CCA-MLE (d = {dc})", S1, s["c_index/cca-mle/multimodal"]["mean"], s["c_index/cca-mle/multimodal"]["std"]))
        items.append((f"Ours: EM (d = {de})", S2, s["c_index/em/multimodal"]["mean"], s["c_index/em/multimodal"]["std"]))
    fig, ax = plt.subplots(figsize=(6.4, 3.3))
    y = np.arange(len(items))[::-1]
    for yi, (lab, col, m, sd) in zip(y, items):
        ax.barh(yi, m, 0.6, color=col, edgecolor=SURF, lw=1.2, zorder=3)
        ax.errorbar(m, yi, xerr=sd, fmt="none", ecolor=INK2, elinewidth=1.0, capsize=3, zorder=4)
        ax.text(m + sd + 0.6, yi, f"{m:.1f}", va="center", fontsize=8, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items])
    ax.set_xlim(40, 85); ax.set_xlabel("test C-index (%), multimodal / concatenated")
    ax.grid(axis="x", color=MUTED, lw=0.5, alpha=0.5); ax.set_axisbelow(True)
    eng = "authors' code" if pb["engine"] == "paper" else "authors' trained models, our Cox"
    ax.set_title("Paper baselines (rerun) vs our latent-variable PGM, mean ± sd over 5 folds\n"
                 f"paper rows: {eng}; our rows: paper's 50 + 50 PCA input, d chosen on validation",
                 fontsize=9, loc="left", color=INK)
    fig.savefig(os.path.join(out, "fig6_paper_comparison.png"))
    plt.close(fig)


if __name__ == "__main__":
    out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "../../Results")
    j = lambda n: json.load(open(os.path.join(out, n), encoding="utf-8")) if os.path.exists(os.path.join(out, n)) else None
    res = j("m2_results_main.json")
    fig1_graphical_model(out)
    if res:
        fig2_em_convergence(res, out)
        fig3_cindex(res, out)
        fig5_entropy(res, out)
    pp = j("m2_results_paperpca.json")
    if pp:
        fig3_cindex(pp, out, "fig3b_cindex_paper_input.png")
    for tag, name in (("main", "fig4_latent_dim_sweep.png"), ("paperpca", "fig4b_latent_dim_sweep_paper_input.png")):
        sw = j(f"m2_latent_dim_sweep_{tag}.json") or (j("m2_latent_dim_sweep.json") if tag == "main" else None)
        if sw:
            fig4_latent_dim(sw, out, name)
    pb = j("m2_paper_baseline.json")
    if pb:
        fig6_paper_comparison(pb, pp, out)
    print("figures written to", out)
