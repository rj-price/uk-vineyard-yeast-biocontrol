#!/usr/bin/env python
"""
28_figS3_assembly.py -> figures/manuscript/FigS3_assembly_rescue.{svg,png}
Supplementary figure: the highly heterozygous ascomycete genomes (Hanseniaspora, Pichia)
fragment under standard short-read (SPAdes) assembly and are rescued by haplotype-aware
Platanus-allee assembly. (A) BUSCO completeness before/after; (B) contig N50 before/after.
Data: analysis/data/platanus_comparison.csv (26 ascomycete isolates).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs

OUT = "analysis/figures/manuscript"
df = pd.read_csv("analysis/data/platanus_comparison.csv")
df["genus"] = df.species.str.split().str[0]
df = df.sort_values(["genus", "new_busco"]).reset_index(drop=True)

COL = {"Hanseniaspora": "#CC3311", "Pichia": "#0077BB"}
gc = df.genus.map(COL).values
y = np.arange(len(df))


def dumbbell(ax, old, new, xlabel, logx=False):
    for i in range(len(df)):
        ax.plot([old[i], new[i]], [y[i], y[i]], color="#cccccc", lw=1.2, zorder=1)
    ax.scatter(old, y, s=14, color="#999999", zorder=2, label="SPAdes")
    ax.scatter(new, y, s=16, color=gc, zorder=3, label="Platanus-allee")
    ax.set_yticks(y)
    ax.set_yticklabels(df.isolate, fontsize=5.2)
    ax.set_xlabel(xlabel, fontsize=7)
    if logx:
        ax.set_xscale("log")
    ax.tick_params(labelsize=6, length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.margins(y=0.01)


def main():
    fs.apply()
    fig = plt.figure(figsize=(7.4, 6.2))
    axA = fig.add_axes([0.10, 0.09, 0.37, 0.79])
    axB = fig.add_axes([0.60, 0.09, 0.34, 0.79])

    dumbbell(axA, df.old_busco.values, df.new_busco.values, "BUSCO complete (%)")
    axA.set_xlim(0, 100)
    dumbbell(axB, df.old_n50.values, df.new_n50.values, "Contig N50 (bp, log scale)", logx=True)
    axB.set_yticklabels([])

    # legends: assembler (grey vs coloured) + genus colours
    from matplotlib.lines import Line2D
    asm = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#999999", markersize=6, label="SPAdes"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="#555555", markersize=6, label="Platanus-allee")]
    gen = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=6,
                  label=fs.patho_it(g)) for g, c in COL.items()]
    axA.legend(handles=asm, fontsize=6, loc="lower right", frameon=False, title="Assembler",
               title_fontsize=6)
    axB.legend(handles=gen, fontsize=6, loc="lower right", frameon=False,
               title="Platanus-allee point", title_fontsize=6)

    fs.panel_letter(axA, "A", x=-0.30, y=1.01)
    fs.panel_letter(axB, "B", x=-0.06, y=1.01)
    fig.suptitle("Haplotype-aware rescue of heterozygous ascomycete genomes",
                 fontsize=9, y=0.975, x=0.52)
    fig.text(0.52, 0.945, "SPAdes → Platanus-allee, 26 Hanseniaspora & Pichia isolates",
             ha="center", fontsize=6.5, color=fs.MUTED)
    fs.save(fig, f"{OUT}/FigS3_assembly_rescue")
    print("wrote FigS3_assembly_rescue")
    print(f"BUSCO: SPAdes {df.old_busco.mean():.1f}% -> Platanus {df.new_busco.mean():.1f}% "
          f"(mean gain {df.busco_gain.mean():.1f}); N50 {df.old_n50.mean():.0f} -> {df.new_n50.mean():.0f} bp")


if __name__ == "__main__":
    main()
