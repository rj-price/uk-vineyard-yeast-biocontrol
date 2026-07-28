#!/usr/bin/env python
"""
24_fig6_invivo.py  ->  figures/manuscript/Fig6_labtofruit.{svg,png}
Manuscript Figure 6 — from lab to fruit: the in vitro -> in vivo gap.
(A) exploratory detached grape-berry assay, one representative berry per condition
    (10 dpi). Water- and yeast-only berries stayed infection-free (yeast is not itself
    pathogenic); both B. cinerea groups developed visible grey mould. n = 1 per group,
    qualitative (no quantitative disease scoring) -> exploratory, hence the honest gap.
(B) genome-characterised candidate-prioritisation dashboard: top isolates by a composite
    mechanism score, showing pathogen suppression (%MI) and mechanism profile, coloured
    by species. Basis for the rational, mechanism-informed selection the grant proposes.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from PIL import Image
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs
from _palette import spc, it

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
DATA = f"{PROJ}/analysis/data"
IMG = f"{DATA}/invivo"
OUT = f"{PROJ}/analysis/figures/manuscript"

# condition -> (file, display title, infection status, chip colour)
CONDS = [
    ("Water_only.png", "Water only", "no infection", "#2CA25F"),
    ("Yeast_only.png", "Yeast only", "no infection", "#2CA25F"),
    ("Bc_only.png", "$\\it{B.\\ cinerea}$ only", "grey mould", "#8a6d3b"),
    ("Bc_and_yeast.png", "$\\it{B.\\ cinerea}$ + yeast", "grey mould", "#8a6d3b"),
]


def load_photo(path):
    """Load a berry photo as-is (on its original black background)."""
    return np.asarray(Image.open(path).convert("RGB"))


def panel_A(axes, brack_ax):
    for ax, (fn, title, status, ccol) in zip(axes, CONDS):
        ax.imshow(load_photo(f"{IMG}/{fn}"))
        ax.set_title(title, fontsize=7.5, pad=3)
        ax.set_facecolor("black")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color("#333333"); sp.set_linewidth(0.6)
        ax.set_frame_on(True)
        # status chip beneath the image
        ax.text(0.5, -0.08, status, transform=ax.transAxes, ha="center", va="top",
                fontsize=6, color="white",
                bbox=dict(boxstyle="round,pad=0.28", fc=ccol, ec="none"))
    # grouping brackets across the two control panels and the two pathogen panels
    brack_ax.axis("off")
    brack_ax.set_xlim(0, 1); brack_ax.set_ylim(0, 1)
    for x0, x1, lab in [(0.005, 0.495, "Controls (no pathogen)"),
                        (0.505, 0.995, "+ $\\it{B.\\ cinerea}$  (10 dpi)")]:
        brack_ax.plot([x0, x0, x1, x1], [0.35, 0.7, 0.7, 0.35], color=fs.MUTED, lw=0.8)
        brack_ax.text((x0 + x1) / 2, 0.78, lab, ha="center", va="bottom", fontsize=7,
                      color=fs.INK)


PATH_COLS = [("MI_Bcinerea", "$\\it{B.c.}$"), ("MI_Foxysporum", "$\\it{F.o.}$"),
             ("MI_Fgraminearum", "$\\it{F.g.}$")]
MECH_COLS = [("VOC_mean", "VOC"), ("Fe_comp_mean", "Fe comp."),
             ("enzyme_breadth", "enzymes")]


def panel_B(ax, cax):
    d = pd.read_csv(f"{DATA}/candidate_ranking.csv")
    d = d[d.species.notna()].sort_values("score", ascending=False).head(6).reset_index(drop=True)
    cols = PATH_COLS + MECH_COLS
    keys = [c[0] for c in cols]
    M = d[keys].astype(float).values
    # per-column min-max normalisation (mixed units) -> relative strength profile
    Z = np.zeros_like(M)
    for j in range(M.shape[1]):
        col = M[:, j]
        rng = np.nanmax(col) - np.nanmin(col)
        Z[:, j] = 0.0 if rng == 0 else (col - np.nanmin(col)) / rng
    im = ax.imshow(Z, aspect="auto", cmap="Purples", vmin=0, vmax=1)
    im.cmap.set_bad("#f0f0f0")

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], fontsize=6.5)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d.isolate, fontsize=7)
    for yi in range(len(d)):
        for xi in range(len(cols)):
            v = M[yi, xi]
            if np.isnan(v):
                continue
            txt = f"{v:.0f}" if keys[xi] != "enzyme_breadth" else f"{int(v)}"
            ax.text(xi, yi, txt, ha="center", va="center", fontsize=6,
                    color="white" if Z[yi, xi] > 0.6 else fs.INK)
        # species colour chip + score to the left / right
        ax.add_patch(Rectangle((-0.135, yi - 0.42), 0.03, 0.84, color=spc(d.species[yi]),
                     transform=ax.get_yaxis_transform(), clip_on=False, lw=0))
    # separator between pathogen-suppression block and mechanism block
    ax.axvline(len(PATH_COLS) - 0.5, color="white", lw=2.4)
    ax.text((len(PATH_COLS) - 1) / 2, -0.9, "pathogen %MI", ha="center", fontsize=6.5,
            color=fs.MUTED, style="italic")
    ax.text(len(PATH_COLS) + (len(MECH_COLS) - 1) / 2, -0.9, "mechanism profile",
            ha="center", fontsize=6.5, color=fs.MUTED, style="italic")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    # score column to the right of the heatmap
    x_sc = len(cols) - 0.5 + 0.5
    for yi in range(len(d)):
        ax.text(len(cols) - 0.3, yi, f"{d.score[yi]:.1f}", va="center", ha="left",
                fontsize=6.5, fontweight="bold", color=fs.INK,
                transform=ax.transData, clip_on=False)
    ax.text(len(cols) - 0.3, -0.9, "score", ha="left", fontsize=6.5, color=fs.MUTED,
            style="italic")

    cb = plt.colorbar(im, cax=cax, ticks=[0, 1])
    cb.set_label("relative\n(per column)", fontsize=5.5)
    cb.ax.set_yticklabels(["min", "max"]); cb.ax.tick_params(labelsize=5.5, length=2)

    # species legend (only the species present among the top candidates)
    from matplotlib.lines import Line2D
    seen = list(dict.fromkeys(d.species))
    handles = [Line2D([0], [0], marker="s", ls="", ms=6, color=spc(s), label=it(s))
               for s in seen]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.14),
              frameon=False, fontsize=6, ncol=2, handletextpad=0.3, columnspacing=1.0)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    fig = plt.figure(figsize=(7.2, 6.6))

    # Panel A: four berry photos in a row (top), with a grouping-bracket strip above
    brack = fig.add_axes([0.06, 0.905, 0.90, 0.05])
    imgw, gap, x0 = 0.205, 0.03, 0.06
    axesA = []
    for i in range(4):
        axesA.append(fig.add_axes([x0 + i * (imgw + gap), 0.60, imgw, 0.29]))
    panel_A(axesA, brack)

    # Panel B: prioritisation dashboard (bottom)
    axB = fig.add_axes([0.16, 0.14, 0.60, 0.30])
    cax = fig.add_axes([0.86, 0.20, 0.018, 0.18])
    panel_B(axB, cax)

    fig.text(0.02, 0.955, "A", fontsize=12, fontweight="bold", va="top")
    fig.text(0.02, 0.47, "B", fontsize=12, fontweight="bold", va="top")

    fig.text(0.06, 0.51,
             "Exploratory assay — one representative berry per condition, "
             "no quantitative disease scoring.\nStrong in vitro antagonism did not "
             "translate to berry protection: the in vitro→in vivo gap this project "
             "is designed to close.",
             fontsize=6.2, color=fs.MUTED, va="top")

    fs.save(fig, f"{OUT}/Fig6_labtofruit")
    print("wrote Fig6_labtofruit")


if __name__ == "__main__":
    main()
