#!/usr/bin/env python
"""
22_fig4_genomic.py  ->  figures/manuscript/Fig4_genomic.{svg,png}
Manuscript Figure 4 — genomic basis of biocontrol traits (7 species-representative
mined genomes; the 5 conspecific comparative genomes are used in Fig 5).
(A) cell-wall-degrading CAZyme inventory (counts per class).
(B) protease arsenal: total vs SignalP-secreted (the secretome refines it).
(C) antiSMASH secondary-metabolite clusters incl. pulcherrimin CDPS (Metschnikowia only).
(D) KT395 killer toxin: present in all P. kluyveri, absent elsewhere.
"""
import os, sys
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
            "Metschnikowia pulcherrima": "Metschnikowia cf. pulcherrima",
            "Metschnikowia ziziphicola": "Metschnikowia cf. pulcherrima",
            "Pichia sp": "Pichia kluyveri", "Pichia sp.": "Pichia kluyveri",
            "Rhodotorula babjevae": "Rhodotorula cf. babjevae",
            "Sporobolomyces roseus": "Sporobolomyces roseus",
            "Vishniacozyma carnescens": "Vishniacozyma carnescens"}.get(raw, raw)


def order(inv):
    inv["sp"] = [genome_species(i, s) for i, s in zip(inv.isolate, inv.species)]
    rank = {s: i for i, s in enumerate(SP_ORDER)}
    inv["_r"] = inv.sp.map(lambda s: rank.get(s, 99))
    return inv.sort_values(["_r", "isolate"]).reset_index(drop=True)


def chips(ax, ys, sps, x=-0.30, w=0.05):
    for y, sp in zip(ys, sps):
        ax.add_patch(Rectangle((x, y - 0.42), w, 0.84, color=spc(sp),
                     transform=ax.get_yaxis_transform(), clip_on=False, lw=0))


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
    ax.set_title("Cell-wall-degrading enzymes (CAZymes)", fontsize=8, pad=4)
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
    ax.set_xlabel("proteases (MEROPS)", fontsize=7)
    ax.tick_params(axis="y", length=0, pad=8)
    chips(ax, y, inv.sp.values, x=-0.16, w=0.03)
    ax.set_title("Protease arsenal: secreted of total", fontsize=8, pad=4)
    h = [Rectangle((0, 0), 1, 1, color="#777", alpha=0.28),
         Rectangle((0, 0), 1, 1, color="#777")]
    ax.legend(h, ["total", "secreted (SignalP)"], loc="upper right", frameon=False,
              fontsize=6, handlelength=1.1, handletextpad=0.4)


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
    ax.set_xlabel("antiSMASH clusters", fontsize=7)
    ax.tick_params(axis="y", length=0, pad=8)
    chips(ax, y, inv.sp.values, x=-0.16, w=0.03)
    ax.set_title("Secondary-metabolite clusters", fontsize=8, pad=4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4, frameon=False,
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
            "absent in all 50 other isolates (only GH18-domain cross-hits)",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.6,
            color=fs.MUTED, style="italic")


SP_DISPLAY = {"Hanseniaspora uvarum": "Hanseniaspora uvarum",
              "Metschnikowia sp.": "Metschnikowia cf. pulcherrima",
              "Pichia sp.": "Pichia kluyveri",
              "Rhodotorula babjevae": "Rhodotorula cf. babjevae",
              "Sporobolomyces roseus": "Sporobolomyces roseus",
              "Vishniacozyma carnescens": "Vishniacozyma carnescens"}


def species_legend(ax, inv):
    present = [s for s in SP_ORDER if s in set(inv.sp)]
    handles = [Line2D([0], [0], marker="s", ls="", ms=7, color=spc(s),
                      label=it(SP_DISPLAY.get(s, s))) for s in present]
    ax.axis("off")
    ax.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=7,
              handletextpad=0.3, columnspacing=1.6, labelspacing=0.6)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    inv = order(pd.read_csv(f"{DATA}/mining_gene_inventory.csv"))
    sec = pd.read_csv(f"{DATA}/secretome_inventory.csv")

    fig = plt.figure(figsize=(7.4, 7.6))
    axA = fig.add_axes([0.11, 0.58, 0.34, 0.34])
    axB = fig.add_axes([0.60, 0.58, 0.34, 0.34])
    axC = fig.add_axes([0.11, 0.10, 0.34, 0.30])
    axD = fig.add_axes([0.62, 0.14, 0.30, 0.24])

    axLeg = fig.add_axes([0.08, 0.45, 0.88, 0.07])
    panel_A(axA, inv)
    panel_B(axB, sec, inv)
    panel_C(axC, inv)
    panel_D(axD)
    species_legend(axLeg, inv)
    fs.panel_letter(axA, "A", x=-0.18, y=1.05)
    fs.panel_letter(axB, "B", x=-0.14, y=1.05)
    fs.panel_letter(axC, "C", x=-0.18, y=1.06)
    fs.panel_letter(axD, "D", x=-0.20, y=1.08)

    fs.save(fig, f"{OUT}/Fig4_genomic")
    print("wrote Fig4_genomic; genomes:", list(inv.isolate))


if __name__ == "__main__":
    main()
