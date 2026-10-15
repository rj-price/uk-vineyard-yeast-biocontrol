#!/usr/bin/env python
"""
25_figS1_global.py -> figures/manuscript/Fig5_global.{svg,png}
Manuscript Figure 5 (promoted from former Fig S1) — global genomic context of the
two lead biocontrol lineages.
(A) ANI within the Rhodotorula babjevae/graminis complex: the Kent isolates form a
    distinct third species-level cluster (incl. misannotated genomes from China &
    Australia), separate from both the R. babjevae and R. graminis type strains.
(B) BUSCO-supermatrix phylogeny (midpoint-rooted phylogram) of the 22 Kent H. uvarum
    and the global H. uvarum genomes (re-deposits of the same strain removed, see
    global_context/meta/duplicate_assemblies.tsv); Kent tips filled if the isolate
    suppressed B. cinerea in the dual-culture screen.
"""
import os, sys, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from Bio import Phylo
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
GC = f"{PROJ}/analysis/global_context"
OUT = f"{PROJ}/analysis/figures/manuscript"

# ---------------- Panel A: Rhodotorula complex ANI ----------------
# curated genomes (label, accession-or-isolate, cluster)
RHO_ROWS = [
    ("UY1  (Kent, UK)", "UY1", "lineage"),
    ("TY27  (Kent, UK)", "TY27", "lineage"),
    ("UY3  (Kent, UK)", "UY3", "lineage"),
    ("JJ10.1  (China)", "GCA_026119225.1", "lineage"),
    ("DBVPG 8058", "GCA_920103735.2", "lineage"),
    ("MAG  (Australia)", "GCA_937859235.1", "lineage"),
    ("CBS 7808$^{T}$  babjevae", "GCA_920104245.2", "babjevae"),
    ("JCM9279  babjevae", "GCA_050403105.1", "babjevae"),
    ("WP1  graminis", "GCA_001329695.1", "graminis"),
    ("JCM3775  graminis", "GCA_050403165.1", "graminis"),
    ("CBS 20  glutinis", "GCA_920103745.3", "glutinis"),
]
# cluster colours kept clear of the yellow-to-red ANI scale; the Kent lineage uses the
# Rhodotorula colour from Figs 1-4
CLUST_COL = {"lineage": "#009E73", "babjevae": "#0077BB", "graminis": "#AA3377",
             "glutinis": "#999999"}


def load_ani(tag):
    names = pd.read_csv(f"{GC}/{tag}_names.csv")
    p2n = dict(zip(names.path.map(os.path.abspath), names.name))
    a = pd.read_csv(f"{GC}/{tag}_ani.tsv", sep="\t", header=None,
                    names=["q", "r", "ani", "frag", "tot"])
    a["qn"] = a.q.map(lambda p: p2n.get(os.path.abspath(p)))
    a["rn"] = a.r.map(lambda p: p2n.get(os.path.abspath(p)))
    return a


def panel_A(ax, cax):
    a = load_ani("rhodo")
    ids = [r[1] for r in RHO_ROWS]
    labs = [r[0] for r in RHO_ROWS]
    piv = (a[a.qn.isin(ids) & a.rn.isin(ids)]
           .pivot_table(index="qn", columns="rn", values="ani")
           .reindex(index=ids, columns=ids))
    M = piv.values.astype(float)
    im = ax.imshow(M, cmap="YlOrRd", vmin=88, vmax=100, aspect="equal")
    im.cmap.set_bad("#eeeeee")
    ax.set_xticks(range(len(ids))); ax.set_yticks(range(len(ids)))
    ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=5.6)
    ax.set_yticklabels(labs, fontsize=5.6)
    for i in range(len(ids)):
        for j in range(len(ids)):
            v = M[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=4.6,
                        color="white" if v >= 96 else fs.INK)
    # colour the tick labels by cluster + draw cluster boxes
    clusters = [r[2] for r in RHO_ROWS]
    for tk, cl in zip(ax.get_yticklabels(), clusters):
        tk.set_color(CLUST_COL[cl])
    for tk, cl in zip(ax.get_xticklabels(), clusters):
        tk.set_color(CLUST_COL[cl])
    # boxes around contiguous clusters
    i = 0
    while i < len(clusters):
        j = i
        while j + 1 < len(clusters) and clusters[j + 1] == clusters[i]:
            j += 1
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), j - i + 1, j - i + 1, fill=False,
                     edgecolor=CLUST_COL[clusters[i]], lw=1.8))
        i = j + 1
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("$\\it{Rhodotorula\\ babjevae/graminis}$ complex",
                 fontsize=7.5, pad=4)
    cb = plt.colorbar(im, cax=cax, ticks=[88, 95, 100], extend="min")
    cb.set_label("ANI (%)", fontsize=6); cb.ax.tick_params(labelsize=5.5, length=2)
    cb.ax.axhline(95, color="#111", lw=1.0)


def rho_legend_handles():
    return [Line2D([0], [0], marker="s", ls="", ms=6, color=CLUST_COL[k],
            label={"lineage": "Kent lineage (+ China, Australia)",
                   "babjevae": "$\\it{R.\\ babjevae}$ type",
                   "graminis": "$\\it{R.\\ graminis}$",
                   "glutinis": "$\\it{R.\\ glutinis}$ (outgroup)"}[k])
            for k in ["lineage", "babjevae", "graminis", "glutinis"]]


# ---------------- Panel B: H. uvarum global phylogeny ----------------
COUNTRY_COL = {
    "Australia": "#EE7733", "China": "#117733", "Italy": "#009988",
    "USA": "#0077BB", "South Korea": "#AA3377", "Belgium": "#CCBB44",
    "Japan": "#EE3377", "Ukraine": "#33BBEE", "unknown": "#BBBBBB",
}


def layout(tree):
    """Phylogram coordinates: x = root-to-node distance, y = leaf order."""
    leaves = tree.get_terminals()
    yy = {lf: i for i, lf in enumerate(leaves)}
    xx = {}

    def dist(cl, d):
        xx[cl] = d
        for c in cl.clades:
            dist(c, d + (c.branch_length or 0))
    dist(tree.root, 0.0)

    def yp(cl):
        if cl.is_terminal():
            return yy[cl]
        y = np.mean([yp(c) for c in cl.clades]); yy[cl] = y; return y
    yp(tree.root)
    return xx, yy, leaves


def panel_B(ax):
    tree = Phylo.read(f"{GC}/huvarum_phylo/supermatrix/huvarum_tree.contree", "newick")
    dup = set(pd.read_csv(f"{GC}/meta/duplicate_assemblies.tsv", sep="\t").acc)
    for lf in [t for t in tree.get_terminals() if t.name in dup]:
        tree.prune(lf)
    tree.root_at_midpoint(); tree.ladderize()
    pre = pd.read_csv(f"{PROJ}/analysis/data/prelim_screen.csv").set_index("isolate")
    bc_active = set(pre.index[pre.GI_Bc_ave.fillna(0) > 0])
    names = pd.read_csv(f"{GC}/huvarum_names.csv")
    kent = set(names[names.kind == "isolate"].name)
    annot = pd.read_csv(f"{GC}/meta/huvarum_annot.csv")
    cc = dict(zip(annot.acc, annot.country)); src = dict(zip(annot.acc, annot.source))
    xx, yy, leaves = layout(tree)
    xmax = max(xx.values()); dx = xmax * 0.03

    for cl in tree.get_nonterminals():
        cys = [yy[c] for c in cl.clades]
        ax.plot([xx[cl], xx[cl]], [min(cys), max(cys)], color="#666", lw=0.6, zorder=1)
        for c in cl.clades:
            ax.plot([xx[cl], xx[c]], [yy[c], yy[c]], color="#666", lw=0.6, zorder=1)
        if cl.confidence is not None and cl is not tree.root and float(cl.confidence) >= 95:
            ax.scatter(xx[cl], yy[cl], s=4, color="#333", zorder=3, edgecolor="none")

    for lf in leaves:
        if lf.name in kent:
            act = lf.name in bc_active
            ax.scatter(xmax + dx, yy[lf], s=15, color="#CC3311" if act else "white",
                       edgecolor="#CC3311", lw=0.7, zorder=4, clip_on=False)
            ax.text(xmax + 2.5 * dx, yy[lf], f"{lf.name}  (Kent, UK)", va="center", ha="left",
                    fontsize=4.3, color="#CC3311", fontweight="bold", clip_on=False)
        else:
            co = cc.get(lf.name, "unknown")
            grape = "grape" in str(src.get(lf.name, ""))
            ax.scatter(xmax + dx, yy[lf], s=15, marker="D" if grape else "o",
                       color=COUNTRY_COL.get(co, "#BBBBBB"), edgecolor="none",
                       zorder=4, clip_on=False)
            tag = " ◆ grape/wine" if grape else ""
            ax.text(xmax + 2.5 * dx, yy[lf], f"{co}{tag}", va="center", ha="left", fontsize=4.3,
                    color=fs.INK, clip_on=False)

    for lf in leaves:   # dotted guide from tip to the aligned label column
        ax.plot([xx[lf], xmax + dx], [yy[lf], yy[lf]], color="#CCC", lw=0.3, ls=":", zorder=0)
    sb = 0.001   # scale bar, substitutions per site
    ax.plot([0, sb], [len(leaves) + 0.2] * 2, color=fs.INK, lw=0.8, clip_on=False)
    ax.text(sb / 2, len(leaves) + 0.6, f"{sb} substitutions/site", ha="center", va="top",
            fontsize=5, clip_on=False)
    ax.set_xlim(-xmax * 0.02, xmax * 1.45)
    ax.set_ylim(len(leaves) - 0.5, -1.5)
    ax.axis("off")
    ax.set_title("$\\it{Hanseniaspora\\ uvarum}$, Kent and global isolates",
                 fontsize=7.5, pad=2)
    handles = [Line2D([0], [0], marker="o", ls="", ms=5, color="#CC3311",
               label="Kent, B. cinerea antagonist"),
               Line2D([0], [0], marker="o", ls="", ms=5, mfc="white", mec="#CC3311",
               label="Kent, not active in screen")]
    handles += [Line2D([0], [0], marker="D", ls="", ms=5, color="#555",
                label="global, grape/wine"),
                Line2D([0], [0], marker="o", ls="", ms=5, color="#555",
                label="global, other source"),
                Line2D([0], [0], marker="o", ls="", ms=4, color="#333",
                label="UFBoot ≥ 95")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(1.0, 0.0),
              frameon=False, fontsize=5.5, handletextpad=0.3)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    fig = plt.figure(figsize=(7.4, 8.8))
    axA = fig.add_axes([0.18, 0.66, 0.42, 0.29])
    caxA = fig.add_axes([0.63, 0.68, 0.015, 0.21])
    axB = fig.add_axes([0.02, 0.02, 0.55, 0.50])
    panel_A(axA, caxA)
    panel_B(axB)
    fig.legend(handles=rho_legend_handles(), loc="upper left",
               bbox_to_anchor=(0.72, 0.90), frameon=False, fontsize=6,
               handletextpad=0.3, title="Panel A clusters",
               title_fontsize=6.5)
    fig.text(0.02, 0.98, "A", fontsize=12, fontweight="bold", va="top")
    fig.text(0.02, 0.54, "B", fontsize=12, fontweight="bold", va="top")
    fs.save(fig, f"{OUT}/Fig5_global")
    print("wrote Fig5_global (promoted from former Fig S1)")


if __name__ == "__main__":
    main()
