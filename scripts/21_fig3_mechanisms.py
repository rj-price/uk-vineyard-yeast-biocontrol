#!/usr/bin/env python
"""
21_fig3_mechanisms.py  ->  figures/manuscript/Fig3_mechanisms.{svg,png}
Manuscript Figure 3 — complementary mechanistic assays on representative isolates
selected per genus (n = 11).
(A) VOC-mediated %MI (volatile assay), per pathogen.
(B) iron competition: %MI at low (5) vs high (20 mg/L) Fe — a leftward shift with
    added Fe indicates iron-dependent antagonism (competition). Strong for the
    Fusaria in Rhodotorula & Metschnikowia; weak/mixed for B. cinerea.
(C) lytic-enzyme relative activity (REA) heatmap, 7 enzymes.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs
from _palette import spc, it, SP_ORDER

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
DATA = f"{PROJ}/analysis/data"
OUT = f"{PROJ}/analysis/figures/manuscript"

ENZ = ["enz_Amylase", "enz_Beta-glucosidase", "enz_Cellulase", "enz_Chitinase",
       "enz_Esterase", "enz_Lipase", "enz_Protease"]
ENZ_LAB = [e.replace("enz_", "").replace("Beta-", "β-") for e in ENZ]


def norm_species(s):
    if not isinstance(s, str) or not s:
        return "Hanseniaspora uvarum"   # RY6: sequencing failed, Hanseniaspora by ITS
    if s.startswith("Hanseniaspora"): return "Hanseniaspora uvarum"
    if s.startswith("Pichia"): return "Pichia sp."
    if s.startswith("Metschnikowia"): return "Metschnikowia sp."
    if "babjevae" in s: return "Rhodotorula babjevae"
    return "Other"


def load():
    m = pd.read_csv(f"{DATA}/phenotype_matrix.csv")
    mm = pd.read_csv(f"{DATA}/isolate_master.csv")[["isolate", "species_genome"]]
    d = m.merge(mm, on="isolate", how="left")
    d = d[d[ENZ].notna().any(axis=1)].copy()
    d["sp"] = d.species_genome.map(norm_species)
    d["nogenome"] = d.species_genome.isna()
    rank = {s: i for i, s in enumerate(SP_ORDER)}
    d["_r"] = d.sp.map(lambda s: rank.get(s, 99))
    d = d.sort_values(["_r", "isolate"]).reset_index(drop=True)
    return d


def species_chips(ax, ys, sps):
    """Species colour chip to the LEFT of each y-tick label (clears the text)."""
    for y, sp in zip(ys, sps):
        ax.add_patch(Rectangle((-0.185, y - 0.42), 0.035, 0.84, color=spc(sp),
                     transform=ax.get_yaxis_transform(), clip_on=False, lw=0))


def panel_A(ax, d):
    P = [("B. cinerea", "VOC_Bcinerea"), ("F. oxysporum", "VOC_Foxysporum"),
         ("F. graminearum", "VOC_Fgraminearum")]
    n = len(d); h = 0.26
    ys = np.arange(n)
    for k, (pth, col) in enumerate(P):
        vals = d[col].values
        off = (k - 1) * h
        ax.barh(ys + off, np.nan_to_num(vals), height=h, color=fs.PATHO_COL[pth],
                zorder=3, label=fs.patho_it(pth))
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r.isolate}*" if r.nogenome else r.isolate
                        for _, r in d.iterrows()], fontsize=6)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlabel("VOC-mediated %MI", fontsize=7)
    ax.set_title("Volatile antagonism", fontsize=8)
    ax.axvline(0, color="#999", lw=0.6)
    ax.legend(loc="lower right", frameon=False, fontsize=6, handlelength=1,
              handletextpad=0.4)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=10)
    species_chips(ax, np.arange(n), d.sp.values)


def panel_B(ax, d):
    """iron competition dumbbells: low (5) filled -> high (20 mg/L) open.
    Leftward shift (added Fe relieves inhibition) = iron-dependent antagonism."""
    pairs = [("B. cinerea", "Fe_Bcinerea_5mgLFe", "Fe_Bcinerea_20mgLFe"),
             ("F. oxysporum", "Fe_Foxysporum_5mgLFe", "Fe_Foxysporum_20mgLFe"),
             ("F. graminearum", "Fe_Fgraminearum_5mgLFe", "Fe_Fgraminearum_20mgLFe")]
    rows = []
    for pth, c5, c20 in pairs:
        for _, r in d.iterrows():
            if not np.isnan(r[c5]) and not np.isnan(r[c20]):
                rows.append((pth, r.isolate, r.sp, r.nogenome, r[c5], r[c20]))
    # order: group by pathogen, then species order within
    rank = {s: i for i, s in enumerate(SP_ORDER)}
    porder = {"B. cinerea": 0, "F. oxysporum": 1, "F. graminearum": 2}
    rows.sort(key=lambda t: (porder[t[0]], rank.get(t[2], 99), t[1]))
    y = 0; yticks = []; ylabs = []; group_y = {}
    for pth, iso, sp, ng, v5, v20 in rows:
        col = fs.PATHO_COL[pth]
        ax.plot([v5, v20], [y, y], color=col, lw=1.2, zorder=2, alpha=0.7)
        ax.scatter(v5, y, s=26, color=col, zorder=3, edgecolor="none")          # low Fe
        ax.scatter(v20, y, s=26, facecolor="white", edgecolor=col, lw=1.1, zorder=3)  # high Fe
        yticks.append(y); ylabs.append(f"{iso}*" if ng else iso)
        group_y.setdefault(pth, []).append(y)
        y += 1
    # pathogen group labels on the right
    for pth, yy in group_y.items():
        ax.text(ax.get_xlim()[1] if False else 62, np.mean(yy), fs.patho_it(pth),
                color=fs.PATHO_COL[pth], fontsize=6.5, va="center", ha="left",
                rotation=0)
        ax.plot([60.5, 60.5], [min(yy) - 0.3, max(yy) + 0.3], color=fs.PATHO_COL[pth],
                lw=2, clip_on=False)
    # direction cue: added Fe relieves Fusarium inhibition (iron competition)
    fy = min(group_y["F. oxysporum"]) - 0.7
    ax.annotate("", xy=(14, fy), xytext=(34, fy),
                arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))
    ax.text(24, fy - 0.5, "+ Fe relieves inhibition\n(iron competition)", fontsize=5.5,
            ha="center", va="bottom", color="#555", style="italic")
    ax.set_yticks(yticks); ax.set_yticklabels(ylabs, fontsize=6)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(0, 66)
    ax.set_xlabel("%MI", fontsize=7)
    ax.set_title("Iron competition (± Fe)", fontsize=8)
    ax.tick_params(axis="y", length=0, pad=10)
    species_chips(ax, yticks, [r[2] for r in rows])
    for _, yy in group_y.items():
        pass
    # legend for Fe level
    h = [Line2D([0], [0], marker="o", ls="", ms=6, color="#555", label="5 mg/L Fe (low)"),
         Line2D([0], [0], marker="o", ls="", ms=6, mfc="white", mec="#555",
                label="20 mg/L Fe (high)")]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=6,
              handletextpad=0.3, bbox_to_anchor=(1.0, -0.02))


def panel_C(ax, d):
    M = d[ENZ].astype(float).values
    im = ax.imshow(M, aspect="auto", cmap="Purples", vmin=1, vmax=3.7)
    ax.set_xticks(range(len(ENZ)))
    ax.set_xticklabels(ENZ_LAB, rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f"{r.isolate}*" if r.nogenome else r.isolate
                        for _, r in d.iterrows()], fontsize=6)
    for yi in range(len(d)):
        for xi in range(len(ENZ)):
            v = M[yi, xi]
            ax.text(xi, yi, f"{v:.1f}", ha="center", va="center", fontsize=5,
                    color="white" if v > 2.6 else fs.INK)
    ax.tick_params(length=0, pad=10)
    species_chips(ax, np.arange(len(d)), d.sp.values)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Lytic-enzyme activity (REA)", fontsize=8, pad=4)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03, ticks=[1, 2, 3])
    cb.set_label("relative activity", fontsize=6)
    cb.ax.tick_params(labelsize=6, length=2)


def species_legend(ax, d):
    present = [s for s in SP_ORDER if s in set(d.sp)]
    handles = [Line2D([0], [0], marker="s", ls="", ms=6, color=spc(s), label=it(s))
               for s in present]
    handles.append(Line2D([0], [0], marker="", ls="", label="* sequencing failed"))
    ax.axis("off")
    ax.legend(handles=handles, loc="center left", ncol=1, frameon=False, fontsize=6.5,
              handletextpad=0.3, labelspacing=0.7)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    d = load()
    fig = plt.figure(figsize=(7.2, 7.6))
    axA = fig.add_axes([0.10, 0.60, 0.34, 0.32])
    axB = fig.add_axes([0.60, 0.42, 0.34, 0.50])
    axC = fig.add_axes([0.10, 0.08, 0.42, 0.30])
    axL = fig.add_axes([0.60, 0.08, 0.36, 0.24])

    panel_A(axA, d)
    panel_B(axB, d)
    panel_C(axC, d)
    species_legend(axL, d)
    fs.panel_letter(axA, "A", x=-0.20, y=1.04)
    fs.panel_letter(axB, "B", x=-0.18, y=1.02)
    fs.panel_letter(axC, "C", x=-0.16, y=1.06)

    fs.save(fig, f"{OUT}/Fig3_mechanisms")
    print("wrote Fig3_mechanisms; n =", len(d))


if __name__ == "__main__":
    main()
