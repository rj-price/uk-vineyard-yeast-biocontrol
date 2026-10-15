#!/usr/bin/env python
"""
23_fig5_within.py  ->  figures/manuscript/FigS1_clonal_ani.{svg,png}
Supplementary figure (demoted from former Fig 5C) — within-group ANI of the two
conspecific comparison sets, showing the near-clonal UY1/TY27 pair (99.99% ANI)
that underpins the natural-experiment argument. The within-genus phenotype and
matched-gene-content panels now live in main Fig 4 (E/F); the clonal-variant
breakdown is reported in the text.
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs
from _palette import spc, it

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
DATA = f"{PROJ}/analysis/data"
ALL = f"{PROJ}/analysis/allelic"
OUT = f"{PROJ}/analysis/figures/manuscript"

GROUPS = [("Rhodotorula aff. babjevae", "Rhodotorula babjevae",
           ["UY1", "TY27", "TY28", "TY29"], "UY1"),
          ("Hanseniaspora uvarum", "Hanseniaspora uvarum",
           ["RY14", "RY7", "RY25"], "RY14")]
PATHO = [("B. cinerea", "MI_Bcinerea"), ("F. oxysporum", "MI_Foxysporum"),
         ("F. graminearum", "MI_Fgraminearum")]
GENE_METRICS = [("secretome", "n_secreted"), ("terpene\nclusters", "terpene"),
                ("secr.\nprotease", "protease_sec"),
                ("secr. β-1,3-\nglucanase", "b-1,3-glucanase|sec"),
                ("secr.\nchitinase", "Chitinase (GH18/19)|sec")]


def flatlist():
    rows = []
    for palkey, _, isos, antag in GROUPS:
        for iso in isos:
            rows.append((iso, palkey, iso == antag))
    return rows


def panel_A(ax):
    ph = pd.read_csv(f"{DATA}/phenotype_matrix.csv").set_index("isolate")
    rows = flatlist()
    y = 0; yt = []; yl = []; gsep = []
    h = 0.26
    for gi, (palkey, _, isos, antag) in enumerate(GROUPS):
        if gi:
            gsep.append(y - 0.5)
        for iso in isos:
            for k, (pth, col) in enumerate(PATHO):
                v = ph.loc[iso, col]
                ax.barh(y + (1 - k) * h, v, height=h, color=fs.PATHO_COL[pth],
                        zorder=3, label=fs.patho_it(pth) if (gi == 0 and iso == isos[0]) else None)
            lab = iso + (" ◄" if iso == antag else "")
            yt.append(y); yl.append(lab)
            y += 1
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=6.5)
    ax.set_ylim(y - 0.5, -0.5)
    ax.set_xlabel("%MI (dual culture)", fontsize=7)
    ax.set_title("Phenotype: Botrytis varies, Fusarium conserved", fontsize=7.5)
    ax.axvline(0, color="#999", lw=0.5)
    for s in gsep:
        ax.axhline(s, color="#bbb", lw=0.6, ls=(0, (3, 2)))
    ax.legend(loc="lower right", frameon=False, fontsize=6, handlelength=1,
              handletextpad=0.4)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for i, (iso, palkey, isantag) in enumerate(rows):
        ax.add_patch(Rectangle((-0.135, i - 0.42), 0.03, 0.84, color=spc(palkey),
                     transform=ax.get_yaxis_transform(), clip_on=False, lw=0))


def panel_B(ax):
    rc = pd.read_csv(f"{DATA}/comparative_Rhodotorula.csv").set_index("isolate")
    hc = pd.read_csv(f"{DATA}/comparative_Hanseniaspora.csv").set_index("isolate")
    comp = pd.concat([rc, hc])
    rows = flatlist()
    isos = [r[0] for r in rows]
    cols = [m[1] for m in GENE_METRICS]
    M = comp.loc[isos, cols].astype(float).values
    # colour by relative deviation from the WITHIN-GENUS column mean, so genuine
    # flatness reads as white and only the (small) real wiggle shows as faint shade
    n_rho = len(GROUPS[0][2])
    Z = np.zeros_like(M)
    for a, b in [(0, n_rho), (n_rho, len(isos))]:
        for j in range(M.shape[1]):
            blk = M[a:b, j]; mu = blk.mean()
            Z[a:b, j] = 0 if mu == 0 else np.abs(blk - mu) / mu
    im = ax.imshow(Z, aspect="auto", cmap="Purples", vmin=0, vmax=0.30)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([m[0] for m in GENE_METRICS], fontsize=6)
    ax.set_yticks(range(len(isos)))
    ax.set_yticklabels([r[0] + (" ◄" if r[2] else "") for r in rows], fontsize=6.5)
    for yi in range(len(isos)):
        for xi in range(len(cols)):
            ax.text(xi, yi, f"{M[yi,xi]:.0f}", ha="center", va="center", fontsize=6,
                    color=fs.INK)
        ax.add_patch(Rectangle((-0.14, yi - 0.42), 0.03, 0.84, color=spc(rows[yi][1]),
                     transform=ax.get_yaxis_transform(), clip_on=False, lw=0))
    ax.axhline(len(GROUPS[0][2]) - 0.5, color="#333", lw=1.2)
    ax.tick_params(length=0, pad=6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Matched gene content is flat", fontsize=7.5, pad=4)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03, ticks=[0, 0.15, 0.30])
    cb.set_label("deviation from\ngenus mean", fontsize=5.5)
    cb.ax.set_yticklabels(["0", "15%", "30%"]); cb.ax.tick_params(labelsize=5.5, length=2)


def ani_matrix(ax, genus, isos, seqmap, title, star=None):
    d = pd.read_csv(f"{ALL}/{genus}_ani_within.tsv", sep="\t", header=None,
                    names=["q", "r", "ani", "m", "t"])
    n = lambda p: re.sub(r"_EK.*", "", os.path.basename(p)).split(".")[0]
    d["q"] = d.q.map(n); d["r"] = d.r.map(n)
    sid = {seqmap[i]: i for i in isos}          # seq_id -> isolate
    d = d[d.q.isin(sid) & d.r.isin(sid)]
    piv = d.pivot_table(index="q", columns="r", values="ani")
    order_sid = [seqmap[i] for i in isos]
    piv = piv.reindex(index=order_sid, columns=order_sid)
    M = piv.values
    im = ax.imshow(M, aspect="equal", cmap="YlOrRd", vmin=98.6, vmax=100)
    ax.set_xticks(range(len(isos))); ax.set_xticklabels(isos, fontsize=6, rotation=0)
    ax.set_yticks(range(len(isos))); ax.set_yticklabels(isos, fontsize=6)
    for i in range(len(isos)):
        for j in range(len(isos)):
            v = M[i, j]
            if not np.isnan(v):
                txt = f"{v:.2f}" if v < 100 else "—"
                ax.text(j, i, txt, ha="center", va="center", fontsize=5.5,
                        color="white" if v >= 99.9 else fs.INK,
                        fontweight="bold" if (star and isos[i] in star and isos[j] in star
                                              and i != j) else "normal")
    if star:
        si = [isos.index(s) for s in star]
        ax.add_patch(Rectangle((min(si) - 0.5, min(si) - 0.5), 2, 2, fill=False,
                     edgecolor="#111", lw=1.6))
    ax.set_title(title, fontsize=7, pad=3)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def panel_D(ax):
    v = pd.read_csv(f"{DATA}/clonal_variants.csv")
    cat = {"non-coding": int(((v.context == "intergenic") | (v.context == "upstream") |
                              v.effect.isin(["intron", "indel_intron_or_utr"])).sum()),
           "synonymous": int((v.effect == "synonymous").sum()),
           "protein-\naltering": int(v.effect.isin(["missense", "nonsense", "frameshift",
                                                    "inframe_indel", "stop_lost"]).sum())}
    labs = list(cat); vals = [cat[k] for k in labs]
    colours = ["#cccccc", "#9ecae1", "#08519c"]
    yb = np.arange(len(labs))
    ax.barh(yb, vals, color=colours, zorder=3, height=0.6)
    for yi, val in zip(yb, vals):
        ax.text(val + 6, yi, str(val), va="center", fontsize=7)
    ax.set_yticks(yb); ax.set_yticklabels(labs, fontsize=6.5)
    ax.set_ylim(len(labs) - 0.5, -0.5)
    ax.set_xlim(0, max(vals) * 1.25)
    ax.set_xlabel("variants (UY1 vs TY27)", fontsize=7)
    ax.set_title("Clonal pair: 602 SNPs + 310 indels (99.99% ANI)", fontsize=7.5, pad=4)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(0.0, -0.42, "Biocontrol genes NOT enriched among protein-altering variants:\n"
            "12.1% of hit genes vs 9.8% of proteome (Fisher p = 0.39).\n"
            "⇒ Botrytis antagonism is allelic/regulatory, not gene-content-encoded.",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.8, color=fs.MUTED)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    master = pd.read_csv(f"{DATA}/isolate_master.csv")
    seqmap = dict(zip(master.isolate, master.seq_id))

    fig = plt.figure(figsize=(7.0, 3.6))
    axC1 = fig.add_axes([0.10, 0.20, 0.34, 0.62])
    axC2 = fig.add_axes([0.56, 0.30, 0.26, 0.42])

    im = ani_matrix(axC1, "Rhodotorula", ["UY1", "TY27", "TY28", "TY29"], seqmap,
                    "$\\it{Rhodotorula}$ aff. $\\it{babjevae}$", star=["UY1", "TY27"])
    ani_matrix(axC2, "Hanseniaspora", ["RY14", "RY7", "RY25"], seqmap,
               "$\\it{Hanseniaspora\\ uvarum}$")
    cb = fig.colorbar(im, ax=axC2, fraction=0.06, pad=0.06, ticks=[99, 99.5, 100])
    cb.set_label("ANI (%)", fontsize=6); cb.ax.tick_params(labelsize=6, length=2)

    fig.text(0.5, 0.95, "Within-group ANI (near-clonal UY1/TY27 pair boxed)",
             fontsize=8.5, fontweight="bold", va="top", ha="center")

    fs.save(fig, f"{OUT}/FigS1_clonal_ani")
    print("wrote FigS1_clonal_ani (demoted from former Fig 5C)")


if __name__ == "__main__":
    main()
