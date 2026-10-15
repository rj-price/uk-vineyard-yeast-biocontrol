#!/usr/bin/env python
"""
22_fig4_genomic.py  ->  figures/manuscript/Fig4_genomic.{svg,png}
Manuscript Figure 4 — genomic basis of biocontrol traits, across the 13 mined
antagonist genomes (spanning the four antagonist genera; the non-antagonist
Sporobolomyces/Vishniacozyma representatives are excluded).
(A) polysaccharide-degrading CAZyme inventory (counts per class).
(B) proteases: total vs SignalP-secreted.
(C) antiSMASH biosynthetic gene clusters incl. pulcherrimin CDPS (Metschnikowia only).
(D) within-genus contrasts: per-pathogen %MI (left) joined to matched gene content
    (right: total CAZymes, total and secreted proteases, BGCs), coloured by deviation
    from the genus mean. The KT395 killer-toxin result is stated in the text (the
    former single-purpose panel was dropped in revision round 3).
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs
from _palette import spc, it, SP_ORDER

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
DATA = f"{PROJ}/analysis/data"
OUT = f"{PROJ}/analysis/figures/manuscript"

# within-genus natural-experiment groups for panels E/F (conspecific comparators):
# Botrytis activity varies, Fusarium activity conserved, gene content flat.
GROUPS = [("Rhodotorula aff. babjevae", ["UY1", "TY27", "TY28", "TY29"], "UY1"),
          ("Hanseniaspora uvarum", ["RY14", "RY7", "RY25"], "RY14")]
PATHO = [("B. cinerea", "MI_Bcinerea"), ("F. oxysporum", "MI_Foxysporum"),
         ("F. graminearum", "MI_Fgraminearum")]
# matched gene-content metrics for panel D: the same quantities as panels A-C
GENE_METRICS = [("CAZymes", "n_cazyme"), ("proteases", "prot_tot"),
                ("secreted\nproteases", "prot_sec"), ("BGCs", "bgc")]

CAZ = ["Chitinase (GH18/19)", "b-1,3-glucanase (GH16/17/55/64/81/128)",
       "b-1,6-glucanase (GH30)", "Cellulase (GH5/6/7/12/45)", "Amylase (GH13/15)",
       "Mannanase (GH76)", "Pectinase (GH28/PL1/PL3/PL4)"]
CAZ_LAB = ["Chitinase", "β-1,3-glucanase", "β-1,6-glucanase", "Cellulase",
           "Amylase", "Mannanase", "Pectinase"]

# SM cluster-type palette (categorical, distinct from species/pathogen schemes)
SM_TYPES = [("Terpene (VOC)", "SM:Terpenes (VOCs)", "#117733"),
            ("NRPS", "SM:NRPS", "#DDCC77"),
            ("CDPS / pulcherrimin", "SM:CDPS / pulcherrimin (iron)", "#CC3311"),
            ("Other", "SM:Other", "#BBBBBB")]


def genome_species(iso, raw):
    if raw.startswith("Hanseniaspora"): return "Hanseniaspora uvarum"
    if raw.startswith("Metschnikowia"): return "Metschnikowia sp."
    if raw.startswith("Pichia"): return "Pichia sp."
    if "babjevae" in raw or "Rhodotorula" in raw: return "Rhodotorula babjevae"
    if raw.startswith("Sporobolomyces"): return "Sporobolomyces roseus"
    if "Vishniacozyma" in raw: return "Vishniacozyma carnescens"
    return "Other"


def gid(raw):
    """Genome-based display name for the representative."""
    return {"Hanseniaspora uvarum": "Hanseniaspora uvarum",
            "Metschnikowia pulcherrima": "Metschnikowia pulcherrima",
            "Metschnikowia ziziphicola": "Metschnikowia pulcherrima",
            "Pichia sp": "Pichia kluyveri", "Pichia sp.": "Pichia kluyveri",
            "Rhodotorula babjevae": "Rhodotorula aff. babjevae",
            "Sporobolomyces roseus": "Sporobolomyces roseus",
            "Vishniacozyma carnescens": "Vishniacozyma carnescens"}.get(raw, raw)


def order(inv):
    inv["sp"] = [genome_species(i, s) for i, s in zip(inv.isolate, inv.species)]
    rank = {s: i for i, s in enumerate(SP_ORDER)}
    inv["_r"] = inv.sp.map(lambda s: rank.get(s, 99))
    return inv.sort_values(["_r", "isolate"]).reset_index(drop=True)


def chips(ax, ys, sps, **_):
    fs.label_chips(ax, ys, [spc(sp) for sp in sps])


def ylabels(inv):
    return [f"{r.isolate}  {it(gid(r.species))}" for _, r in inv.iterrows()]


def panel_A(ax, inv):
    M = inv[CAZ].astype(float).values
    im = ax.imshow(M, aspect="auto", cmap="YlGnBu", vmin=0, vmax=M.max())
    ax.set_xticks(range(len(CAZ))); ax.set_xticklabels(CAZ_LAB, rotation=40, ha="right",
                                                       fontsize=6.5)
    ax.set_yticks(range(len(inv)))
    ax.set_yticklabels([r.isolate for _, r in inv.iterrows()], fontsize=6.5)
    for yi in range(len(inv)):
        for xi in range(len(CAZ)):
            v = M[yi, xi]
            ax.text(xi, yi, f"{v:.0f}", ha="center", va="center", fontsize=5.5,
                    color="white" if v > M.max() * 0.6 else fs.INK)
    ax.tick_params(length=0, pad=8)
    for s in ax.spines.values():
        s.set_visible(False)
    chips(ax, np.arange(len(inv)), inv.sp.values, x=-0.16, w=0.03)
    ax.set_title("Polysaccharide-degrading enzymes (CAZymes)", fontsize=8, pad=4)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("gene count", fontsize=6); cb.ax.tick_params(labelsize=6, length=2)


def panel_B(ax, sec, inv):
    s = sec.set_index("isolate").loc[inv.isolate]
    tot = s["MEROPS_protease | total"].values
    scr = s["MEROPS_protease | secr"].values
    y = np.arange(len(inv))
    for i, (_, r) in enumerate(inv.iterrows()):
        c = spc(r.sp)
        ax.barh(y[i], tot[i], height=0.62, color=c, alpha=0.28, zorder=2)   # total
        ax.barh(y[i], scr[i], height=0.62, color=c, zorder=3)               # secreted
        ax.text(tot[i] + 4, y[i], f"{int(scr[i])}/{int(tot[i])}", va="center",
                fontsize=5.5, color=fs.INK)
    ax.set_yticks(y); ax.set_yticklabels([r.isolate for _, r in inv.iterrows()],
                                         fontsize=6.5)
    ax.set_ylim(len(inv) - 0.5, -0.5)
    ax.set_xlim(0, 265)
    ax.set_xlabel("proteases", fontsize=7)
    ax.tick_params(axis="y", length=0, pad=8)
    chips(ax, y, inv.sp.values, x=-0.16, w=0.03)
    ax.set_title("Total and secreted proteases", fontsize=8, pad=4)
    h = [Rectangle((0, 0), 1, 1, color="#777", alpha=0.28),
         Rectangle((0, 0), 1, 1, color="#777")]
    ax.legend(h, ["total", "secreted"], loc="upper center",
              bbox_to_anchor=(0.5, -0.17), frameon=False, ncol=2,
              fontsize=6, handlelength=1.1, handletextpad=0.4, columnspacing=1.0)


def panel_C(ax, inv):
    y = np.arange(len(inv))
    left = np.zeros(len(inv))
    for lab, col, c in SM_TYPES:
        vals = inv[col].fillna(0).astype(float).values
        ax.barh(y, vals, left=left, height=0.62, color=c, label=lab, zorder=3,
                edgecolor="white", linewidth=0.6)
        left += vals
    ax.set_yticks(y); ax.set_yticklabels([r.isolate for _, r in inv.iterrows()],
                                         fontsize=6.5)
    ax.set_ylim(len(inv) - 0.5, -0.5)
    ax.set_xlabel("Biosynthetic gene clusters", fontsize=7)
    ax.tick_params(axis="y", length=0, pad=8)
    chips(ax, y, inv.sp.values, x=-0.16, w=0.03)
    ax.set_title("Biosynthetic gene clusters", fontsize=8, pad=4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4, frameon=False,
              fontsize=5.8, handlelength=1, handletextpad=0.4, columnspacing=1.0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)


def panel_D(ax):
    kt = pd.read_csv(f"{DATA}/killer_toxin_hits.csv")
    g = kt[(kt.call == "genuine")].sort_values("isolate")
    y = np.arange(len(g))
    ax.barh(y, g.Bc_dual.values, height=0.6, color=spc("Pichia sp."), zorder=3)
    for yi, (_, r) in zip(y, g.iterrows()):
        ax.text(r.Bc_dual + 0.8, yi, f"{r.Bc_dual:.0f}%", va="center", fontsize=6,
                color=fs.INK)
    ax.set_yticks(y); ax.set_yticklabels(g.isolate, fontsize=6.5)
    ax.set_ylim(len(g) - 0.5, -0.5)
    ax.set_xlim(0, 55)
    ax.set_xlabel(fs.patho_it("B. cinerea") + " %MI", fontsize=7)
    ax.set_title("KT395 killer toxin (" + it("Pichia kluyveri") + ")", fontsize=8, pad=4)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.text(0.5, -0.42, "KT395 present in all 4 P. kluyveri (96% id, full-length);\n"
            "absent in all 51 other isolates (only GH18-domain cross-hits)",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.6,
            color=fs.MUTED, style="italic")


def _within_rows():
    rows = []
    for palkey, isos, antag in GROUPS:
        for iso in isos:
            rows.append((iso, palkey, iso == antag))
    return rows


def panel_E(ax):
    """Within-genus phenotype: %MI per pathogen. Botrytis varies in both contrasts;
    F. graminearum also varies across the Rhodotorula pair, so only F. oxysporum is conserved."""
    ph = pd.read_csv(f"{DATA}/phenotype_matrix.csv").set_index("isolate")
    y = 0; yt = []; yl = []; gsep = []
    h = 0.26
    for gi, (palkey, isos, antag) in enumerate(GROUPS):
        if gi:
            gsep.append(y - 0.5)
        for iso in isos:
            for k, (pth, col) in enumerate(PATHO):
                v = ph.loc[iso, col]
                ax.barh(y + (1 - k) * h, v, height=h, color=fs.PATHO_COL[pth], zorder=3,
                        label=fs.patho_it(pth) if (gi == 0 and iso == isos[0]) else None)
            yt.append(y); yl.append(iso)
            y += 1
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=6.5)
    ax.set_ylim(y - 0.5, -0.5)
    ax.set_xlabel("%MI", fontsize=7)
    ax.set_title("Antagonism", fontsize=7, pad=4)
    ax.axvline(0, color="#999", lw=0.5)
    for s in gsep:
        ax.axhline(s, color="#bbb", lw=0.6, ls=(0, (3, 2)))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False,
              fontsize=6, handlelength=1, handletextpad=0.4, ncol=1)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    rows = _within_rows()
    fs.label_chips(ax, yt, [spc(r[1]) for r in rows])


def panel_F(ax, cax):
    """Matched gene content for the panel-E rows (no labels: shares E's rows)."""
    inv = pd.read_csv(f"{DATA}/mining_gene_inventory.csv").set_index("isolate")
    sec = pd.read_csv(f"{DATA}/secretome_inventory.csv").set_index("isolate")
    comp = pd.DataFrame({"n_cazyme": inv["n_cazyme"],
                         "prot_tot": sec["MEROPS_protease | total"],
                         "prot_sec": sec["MEROPS_protease | secr"],
                         "bgc": inv["n_SM_clusters"]})
    rows = _within_rows()
    isos = [r[0] for r in rows]
    cols = [m[1] for m in GENE_METRICS]
    M = comp.loc[isos, cols].astype(float).values
    n_rho = len(GROUPS[0][1])
    Z = np.zeros_like(M)
    for a, b in [(0, n_rho), (n_rho, len(isos))]:
        for j in range(M.shape[1]):
            blk = M[a:b, j]; mu = blk.mean()
            dev = np.abs(blk - mu)
            # a difference of one gene is within annotation noise: left unshaded (as Fig S2)
            Z[a:b, j] = 0 if mu == 0 else np.where(dev <= 1, 0.0, dev / mu)
    im = ax.imshow(Z, aspect="auto", cmap="Purples", vmin=0, vmax=0.30)
    ax.set_ylim(len(isos) - 0.5, -0.5)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([m[0] for m in GENE_METRICS], fontsize=6, rotation=30, ha="right")
    ax.set_yticks([])
    for yi in range(len(isos)):
        for xi in range(len(cols)):
            ax.text(xi, yi, f"{M[yi,xi]:.0f}", ha="center", va="center", fontsize=6,
                    color="white" if Z[yi, xi] > 0.18 else fs.INK)
    ax.axhline(n_rho - 0.5, color="#bbb", lw=0.6, ls=(0, (3, 2)))
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set_title("Gene count", fontsize=7, pad=4)
    cb = plt.colorbar(im, cax=cax, ticks=[0, 0.15, 0.30])
    cb.set_label("deviation from genus mean", fontsize=5.5)
    cb.ax.set_yticklabels(["0", "15%", "30%"]); cb.ax.tick_params(labelsize=5.5, length=2)


SP_DISPLAY = {"Hanseniaspora uvarum": "Hanseniaspora uvarum",
              "Metschnikowia sp.": "Metschnikowia pulcherrima",
              "Pichia sp.": "Pichia kluyveri",
              "Rhodotorula babjevae": "Rhodotorula aff. babjevae",
              "Sporobolomyces roseus": "Sporobolomyces roseus",
              "Vishniacozyma carnescens": "Vishniacozyma carnescens"}


def species_legend(ax, inv):
    present = [s for s in SP_ORDER if s in set(inv.sp)]
    handles = [Line2D([0], [0], marker="s", ls="", ms=7, color=spc(s),
                      label=it(SP_DISPLAY.get(s, s))) for s in present]
    ax.axis("off")
    ax.legend(handles=handles, loc="center", ncol=2, frameon=False, fontsize=7,
              handletextpad=0.3, columnspacing=1.2, labelspacing=0.6,
              title="mined genome species", title_fontsize=7)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    inv = order(pd.read_csv(f"{DATA}/mining_gene_inventory.csv"))
    sec = pd.read_csv(f"{DATA}/secretome_inventory.csv")

    fig = plt.figure(figsize=(7.6, 8.6))
    # Row 1 (13 mined genomes): A CAZymes | B proteases
    axA = fig.add_axes([0.14, 0.62, 0.30, 0.32])
    axB = fig.add_axes([0.62, 0.62, 0.33, 0.32])
    # Row 2: C biosynthetic gene clusters | D within-genus contrasts (E bars + F heatmap)
    TOP2 = 0.44
    axC = fig.add_axes([0.14, 0.10, 0.30, TOP2 - 0.10])
    axE = fig.add_axes([0.60, 0.21, 0.15, TOP2 - 0.21])
    axF = fig.add_axes([0.765, 0.21, 0.16, TOP2 - 0.21])
    cax = fig.add_axes([0.94, 0.21, 0.012, TOP2 - 0.21])
    axLeg = fig.add_axes([0.52, 0.0, 0.46, 0.07])

    panel_A(axA, inv)
    panel_B(axB, sec, inv)
    panel_C(axC, inv)
    panel_E(axE)
    panel_F(axF, cax)
    species_legend(axLeg, inv)
    fig.text(0.765, TOP2 + 0.045, "Within-genus contrasts", ha="center", fontsize=8)

    fs.panel_letter(axA, "A", x=-0.20, y=1.06)
    fs.panel_letter(axB, "B", x=-0.14, y=1.06)
    fig.text(0.02, TOP2 + 0.03, "C", fontsize=12, fontweight="bold", va="bottom", ha="left")
    fig.text(0.51, TOP2 + 0.03, "D", fontsize=12, fontweight="bold", va="bottom", ha="left")

    fs.save(fig, f"{OUT}/Fig4_genomic")
    print("wrote Fig4_genomic; genomes:", list(inv.isolate))


if __name__ == "__main__":
    main()
