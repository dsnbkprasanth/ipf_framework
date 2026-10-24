"""Figure generation mirroring the reference framework's Figs 1-7, for IPF."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold",
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

# colourblind-safe palette
CB = {"up": "#c0392b", "down": "#2471a3", "neutral": "#7f8c8d",
      "fib": "#e67e22", "epi": "#8e44ad", "mye": "#16a085",
      "endo": "#2980b9", "lym": "#c0392b"}
CT_COLOR = {"Fibroblast": CB["fib"], "Epithelial": CB["epi"], "Myeloid": CB["mye"],
            "Endothelial": CB["endo"], "Lymphoid": CB["lym"]}
ROUTE_COLOR = {"TGF-beta": "#c0392b", "PDGF": "#e67e22", "matrix-integrin": "#8e44ad",
               "chemokine": "#16a085", "cytokine": "#2980b9", "WNT": "#d4ac0d"}


def _save(fig, path, tight=True):
    if tight:
        fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return str(path)


# --------------------------------------------------------------------------- #
def fig1_workflow(outdir: Path, demo: bool = False) -> str:
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.axis("off")
    steps = [
        ("Bulk discovery\ncohorts\n(GSE47460,\nGSE150910)", CB["down"]),
        ("Cross-cohort\nreproducible\nsignature", CB["down"]),
        ("Donor-level\npseudobulk\nlocalization\n(GSE135893)", CB["fib"]),
        ("Curated module\nscoring\n(epithelial /\nfibrotic / TGF-b)", CB["epi"]),
        ("Focused\nligand-receptor\nprioritization", CB["mye"]),
        ("External validation\n(GSE134692 bulk,\nGSE136831 snRNA)", CB["endo"]),
        ("LINCS/L1000\nsafety-aware\nperturbational\nprioritization", CB["up"]),
    ]
    n = len(steps)
    x0, w, gap, y, h = 0.2, 1.35, 0.28, 0.5, 1.4
    for i, (txt, col) in enumerate(steps):
        x = x0 + i * (w + gap)
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.12",
                             linewidth=1.4, edgecolor=col, facecolor=col + "22")
        ax.add_patch(box)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=7.6)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((x + w, y + h / 2),
                                         (x + w + gap, y + h / 2),
                                         arrowstyle="-|>", mutation_scale=13,
                                         linewidth=1.3, color="#555"))
    ax.text(x0, y + h + 0.35,
            "Reproducible computational workflow for IPF transcriptomic integration",
            fontsize=11, fontweight="bold")
    ax.text(x0, y - 0.3,
            "Candidate compounds are perturbational hypotheses, not therapeutic evidence."
            + ("   [DEMO DATA]" if demo else ""),
            fontsize=7.5, style="italic", color="#666")
    ax.set_xlim(0, x0 + n * (w + gap))
    ax.set_ylim(0, 2.6)
    return _save(fig, outdir / "Fig1_workflow.png", tight=False)


# --------------------------------------------------------------------------- #
def fig2_crosscohort(outdir: Path, stage1: dict, stage2: dict, stage7: dict) -> str:
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    ids = stage1["cohort_ids"][:2]
    m = stage1["merged"]

    # A: concordance scatter
    ax = axes[0, 0]
    x, y = m[f"log2FC_{ids[0]}"], m[f"log2FC_{ids[1]}"]
    rep = m["replicated_up"] | m["replicated_down"]
    ax.scatter(x[~rep], y[~rep], s=5, c=CB["neutral"], alpha=0.3, label="other")
    ax.scatter(x[m["replicated_up"]], y[m["replicated_up"]], s=8, c=CB["up"], label="replicated up")
    ax.scatter(x[m["replicated_down"]], y[m["replicated_down"]], s=8, c=CB["down"], label="replicated down")
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    r = np.corrcoef(x, y)[0, 1]
    ax.set_xlabel(f"log2FC {ids[0]}"); ax.set_ylabel(f"log2FC {ids[1]}")
    ax.set_title(f"A  Cross-cohort concordance (r={r:.2f})")
    ax.legend(fontsize=6, loc="upper left")

    # B: replicated counts
    ax = axes[0, 1]
    s = stage1["summary"]
    bars = ["same\ndirection", "replicated\nup", "replicated\ndown"]
    vals = [s["n_same_direction"], s["n_replicated_up"], s["n_replicated_down"]]
    ax.bar(bars, vals, color=[CB["neutral"], CB["up"], CB["down"]])
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("genes"); ax.set_title("B  Reproducible IPF signature")

    # C: leave-one-cohort-out
    ax = axes[1, 0]
    loco = stage7["loco"]
    labels = [f"{r.discovery}\n→{r.validation}" for r in loco.itertuples()]
    ax.bar(labels, loco["same_direction_frac"] * 100, color=CB["fib"])
    for i, v in enumerate(loco["same_direction_frac"] * 100):
        ax.text(i, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 100); ax.set_ylabel("same-direction (%)")
    ax.set_title("C  Leave-one-cohort-out sensitivity")

    # D: pseudobulk DEG counts by cell type
    ax = axes[1, 1]
    deg = stage2["deg_counts"]
    cts = list(deg.index)
    xpos = np.arange(len(cts))
    ax.bar(xpos - 0.2, deg["up"], 0.4, color=CB["up"], label="up")
    ax.bar(xpos + 0.2, deg["down"], 0.4, color=CB["down"], label="down")
    ax.set_xticks(xpos); ax.set_xticklabels(cts, rotation=30, ha="right")
    ax.set_ylabel("disease-associated genes")
    ax.set_title("D  Donor-level pseudobulk DEGs")
    ax.legend(fontsize=7)
    return _save(fig, outdir / "Fig2_crosscohort.png")


# --------------------------------------------------------------------------- #
def fig3_modules(outdir: Path, stage3: dict) -> str:
    shift = stage3["module_shift"].astype(float)
    labels = stage3["labels"]
    fig, ax = plt.subplots(figsize=(9, 5.2))
    data = shift.copy()
    data.index = [labels.get(i, i) for i in data.index]
    vmax = np.nanmax(np.abs(data.values))
    im = ax.imshow(data.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(data.shape[1])); ax.set_xticklabels(data.columns, rotation=30, ha="right")
    ax.set_yticks(range(data.shape[0])); ax.set_yticklabels(data.index)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color="white" if abs(v) > vmax * 0.6 else "black")
    fig.colorbar(im, ax=ax, label="mean disease log2FC (module shift)")
    ax.set_title("Fibro-epithelial module shifts by cell type")
    return _save(fig, outdir / "Fig3_modules.png")


# --------------------------------------------------------------------------- #
def fig4_ligand_receptor(outdir: Path, stage4: dict) -> str:
    tab = stage4["axes"].head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    colors = [ROUTE_COLOR.get(r, "#888") for r in tab["route"]]
    ax.barh(tab["axis"], tab["evidence_score"], color=colors)
    ax.set_xlabel("evidence score"); ax.set_title("Focused ligand-receptor prioritization")
    seen = {}
    for r in tab["route"]:
        seen[r] = ROUTE_COLOR.get(r, "#888")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in seen.values()]
    ax.legend(handles, list(seen), fontsize=7, title="route", loc="lower right")
    return _save(fig, outdir / "Fig4_ligand_receptor.png")


# --------------------------------------------------------------------------- #
def fig5_tgfb(outdir: Path, stage3: dict) -> str:
    corr = stage3["module_corr"].astype(float)
    labels = stage3["labels"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    disp = corr.copy()
    disp.index = [labels.get(i, i) for i in disp.index]
    disp.columns = [labels.get(i, i) for i in disp.columns]
    im = ax.imshow(disp.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(disp.shape[1])); ax.set_xticklabels(disp.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(disp.shape[0])); ax.set_yticklabels(disp.index, fontsize=7)
    fig.colorbar(im, ax=ax, label="Spearman r (across cell types)")
    ax.set_title("A  Module co-behavior")

    ax = axes[1]
    if "tgfb_signaling" in corr.index:
        tgfb = corr["tgfb_signaling"].drop("tgfb_signaling").sort_values()
        names = [labels.get(i, i) for i in tgfb.index]
        ax.barh(names, tgfb.values, color=[CB["up"] if v > 0 else CB["down"] for v in tgfb.values])
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("Spearman r with TGF-b coordinate")
    ax.set_title("B  TGF-b coordinate association")
    return _save(fig, outdir / "Fig5_tgfb.png")


# --------------------------------------------------------------------------- #
def fig6_external(outdir: Path, stage3: dict, stage5: dict, stage7: dict) -> str:
    labels = stage3["labels"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))

    # A external bulk module diff
    ax = axes[0]
    eb = stage5["external_bulk"]["IPF_minus_control"].sort_values()
    names = [labels.get(i, i) for i in eb.index]
    ax.barh(names, eb.values, color=[CB["up"] if v > 0 else CB["down"] for v in eb.values])
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("IPF − control module score"); ax.set_title("A  External bulk (GSE134692)")

    # B snRNA external heatmap
    ax = axes[1]
    sn = stage5["snrna_external"].astype(float)
    sn.index = [labels.get(i, i) for i in sn.index]
    vmax = np.nanmax(np.abs(sn.values))
    im = ax.imshow(sn.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(sn.shape[1])); ax.set_xticklabels(sn.columns, rotation=30, ha="right", fontsize=7)
    ax.set_yticks(range(sn.shape[0])); ax.set_yticklabels(sn.index, fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, label="module shift")
    ax.set_title("B  External snRNA (GSE136831)")

    # C robustness null distribution
    ax = axes[2]
    rc = stage7["random_control"]
    ax.hist(rc["null"], bins=30, color=CB["neutral"], alpha=0.7)
    ax.axvline(rc["observed_concordance"], color=CB["up"], lw=2,
               label=f"observed={rc['observed_concordance']:.2f}\nP={rc['empirical_p_str']}")
    ax.set_xlabel("cross-cohort concordance"); ax.set_ylabel("random sets")
    ax.set_title("C  Size-matched random control"); ax.legend(fontsize=7)
    return _save(fig, outdir / "Fig6_external.png")


# --------------------------------------------------------------------------- #
def fig7_perturbation(outdir: Path, stage6: dict) -> str:
    tab = stage6["perturbation"].copy()
    cat_color = {"prioritized": "#1e8449", "supporting": "#7dcea0",
                 "cautionary": "#e67e22", "deprioritized": "#95a5a6"}
    tab = tab.sort_values("total_score")
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = [cat_color.get(c, "#888") for c in tab["category"]]
    ax.barh(tab["compound"], tab["total_score"], color=colors)
    pos = stage6["positive_control_ranks"]
    for cmpd in pos:
        if cmpd in list(tab["compound"]):
            yi = list(tab["compound"]).index(cmpd)
            ax.text(tab["total_score"].iloc[yi], yi, "  ★ approved antifibrotic",
                    va="center", fontsize=7, color="#154360")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("composite translational score")
    ax.set_title("Safety-aware perturbational prioritization")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in cat_color.values()]
    ax.legend(handles, list(cat_color), fontsize=7, loc="lower right", title="category")
    return _save(fig, outdir / "Fig7_perturbation.png")
