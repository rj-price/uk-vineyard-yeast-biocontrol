#!/usr/bin/env python
"""
19_fig1_screen.py  ->  figures/manuscript/Fig1_screen.{svg,png}
Manuscript Figure 1 — the antagonism screen.
(A) preliminary competition screen (Effect / No-effect) of all screened isolates
    x 5 pathogen strains (B. cinerea x3, F. oxysporum, F. graminearum), rows
    grouped by genome species; Verticillium dahliae omitted (project decision).
(B) quantitative %MI (dual-plate) for the isolates taken forward, as a per-pathogen
    heatmap, isolates grouped by species.
(C) pathogen-specificity of the taken-forward set as an UpSet-style intersection.

Source: gideon_mphil/Combined Data.xlsx (preliminary Effect/No-effect + quantitative
GI%) and data/phenotype_matrix.csv (parsed dual-plate %MI).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs
from _palette import spc, it, SP_ORDER

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
XLSX = f"{PROJ}/gideon_mphil/Combined Data.xlsx"
DATA = f"{PROJ}/analysis/data"
OUT = f"{PROJ}/analysis/figures/manuscript"

# preliminary-screen strain columns kept (Vd dropped), grouped by pathogen
STRAINS = [("B. cinerea", ["Bc21010", "Bc21024", "BcCHF"]),
           ("F. oxysporum", ["FoAJ516"]),
           ("F. graminearum", ["Fg321"])]
STRAIN_LABEL = {"Bc21010": "21010", "Bc21024": "21024", "BcCHF": "CHF",
                "FoAJ516": "AJ516", "Fg321": "321"}


def genome_species():
    m = pd.read_csv(f"{DATA}/isolate_master.csv")
    return dict(zip(m.isolate, m.species_genome.fillna("")))


def norm_species(s):
    """Map genome species names onto the fixed _palette SP_ORDER keys."""
    if not s:
        return "Other"
    if s.startswith("Hanseniaspora"): return "Hanseniaspora uvarum"
    if s.startswith("Pichia"): return "Pichia sp."
    if s.startswith("Metschnikowia"): return "Metschnikowia sp."
    if "babjevae" in s: return "Rhodotorula babjevae"
    if s.startswith("Sporobolomyces roseus"): return "Sporobolomyces roseus"
    if s.startswith("Sporobolomyces"): return "Sporobolomyces sp."
    if "carnescens" in s: return "Vishniacozyma carnescens"
    if "victoriae" in s: return "Sporobolomyces sp."  # rare singletons -> group tail
    if s.startswith("Filobasidium"): return "Filobasidium stepposum"
    return "Other"


def load_prelim():
    # cached from Combined Data.xlsx by the funannotate env (data_science lacks openpyxl)
    d = pd.read_csv(f"{DATA}/prelim_screen.csv")
    allstrains = [c for _, cs in STRAINS for c in cs]
    gsp = genome_species()
    d["species"] = d.isolate.map(lambda i: norm_species(gsp.get(i, "")))
    d["n_eff"] = d[allstrains].sum(axis=1)
    return d, allstrains


def order_isolates(d, strains):
    """Sort by species (SP_ORDER), then by number of positive strains desc."""
    rank = {s: i for i, s in enumerate(SP_ORDER)}
    d = d.assign(_r=d.species.map(lambda s: rank.get(s, 99)))
    return d.sort_values(["_r", "n_eff", "isolate"],
                         ascending=[True, False, True]).reset_index(drop=True)


def panel_A(ax, d, strains):
    n = len(d)
    # lay out strains with a 0.8-col gap between pathogen groups so the
    # single-strain Fo/Fg group labels do not collide
    xpos, layout = 0.0, []   # (x, strain, pathogen)
    group_span = {}
    GAP = 1.7
    for pth, cs in STRAINS:
        x0 = xpos
        for c in cs:
            layout.append((xpos, c, pth)); xpos += 1
        group_span[pth] = (x0, xpos - 1)
        xpos += GAP
    xmax = xpos - GAP
    # species colour sidebar
    for y, sp in enumerate(d.species):
        ax.add_patch(Rectangle((-1.5, y - 0.5), 1.0, 1, color=spc(sp), lw=0))
    # cells
    for x, c, pth in layout:
        for y, v in enumerate(d[c].values):
            if v:
                ax.add_patch(Rectangle((x - 0.44, y - 0.44), 0.88, 0.88, color=fs.INK, lw=0))
            else:
                ax.add_patch(Rectangle((x - 0.44, y - 0.44), 0.88, 0.88,
                                       facecolor="#f2f2f2", edgecolor="#e6e6e6", lw=0.3))
    # pathogen group brackets + labels; stagger label heights so the two
    # single-strain Fusarium names (wider than their column) never collide
    label_y = {"B. cinerea": -1.5, "F. oxysporum": -2.35, "F. graminearum": -1.5}
    for pth, (a, b) in group_span.items():
        ax.plot([a - 0.44, b + 0.44], [-1.0, -1.0], color=fs.PATHO_COL[pth],
                lw=2.6, solid_capstyle="butt", clip_on=False)
        ly = label_y[pth]
        if ly < -1.7:   # leader line from bracket down to a lowered label
            ax.plot([(a + b) / 2, (a + b) / 2], [-1.05, ly + 0.15],
                    color=fs.PATHO_COL[pth], lw=0.5, clip_on=False)
        ax.text((a + b) / 2, ly, fs.patho_it(pth), ha="center", va="bottom",
                fontsize=6.5, color=fs.PATHO_COL[pth])
    ax.set_xticks([x for x, _, _ in layout])
    ax.set_xticklabels([STRAIN_LABEL[c] for _, c, _ in layout], fontsize=6,
                       rotation=40, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(d.isolate, fontsize=4.5)
    ax.set_xlim(-1.7, xmax + 0.5)
    ax.set_ylim(n - 0.5, -2.8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"Preliminary screen ({n} isolates)", fontsize=8, pad=14)


def panel_B(ax, order_iso):
    m = pd.read_csv(f"{DATA}/phenotype_matrix.csv").set_index("isolate")
    cols = {"B. cinerea": "MI_Bcinerea", "F. oxysporum": "MI_Foxysporum",
            "F. graminearum": "MI_Fgraminearum"}
    iso = [i for i in order_iso if i in m.index]
    M = m.loc[iso, list(cols.values())].astype(float)
    im = ax.imshow(M.values, aspect="auto", cmap=fs.MI_CMAP, vmin=0, vmax=65)
    ax.set_xticks(range(3))
    ax.set_xticklabels([fs.patho_it(p) for p in cols], rotation=30, ha="right", fontsize=7)
    ax.set_yticks(range(len(iso)))
    ax.set_yticklabels(iso, fontsize=5)
    # annotate values
    for y in range(len(iso)):
        for x in range(3):
            v = M.values[y, x]
            ax.text(x, y, f"{v:.0f}", ha="center", va="center", fontsize=4.5,
                    color="white" if v > 38 else fs.INK)
    ax.set_title("Dual-plate %MI (taken forward)", fontsize=8, pad=6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label("% mycelial inhibition", fontsize=6)
    cb.ax.tick_params(labelsize=6, length=2)
    return iso


def panel_C(ax_bar, ax_mat):
    """Classic UpSet: vertical count bars on top (ax_bar), membership dot-matrix
    below (ax_mat) with pathogen rows labelled on the left."""
    m = pd.read_csv(f"{DATA}/phenotype_matrix.csv").set_index("isolate")
    P = ["B. cinerea", "F. oxysporum", "F. graminearum"]
    col = {"B. cinerea": "MI_Bcinerea", "F. oxysporum": "MI_Foxysporum",
           "F. graminearum": "MI_Fgraminearum"}
    act = {p: set(m.index[m[col[p]] > 0]) for p in P}
    from itertools import combinations
    combos = []
    for k in (1, 2, 3):
        for cmb in combinations(P, k):
            s = set(m.index)
            for p in P:
                s = s & act[p] if p in cmb else s - act[p]
            if s:
                combos.append((set(cmb), len(s)))
    combos.sort(key=lambda t: -t[1])
    x = np.arange(len(combos))

    # ---- bars ----
    ax_bar.bar(x, [c[1] for c in combos], width=0.6, color="#4a4a4a", zorder=3)
    for xi, (_, cnt) in zip(x, combos):
        ax_bar.text(xi, cnt + 0.4, str(cnt), ha="center", va="bottom", fontsize=7)
    XR = len(combos) + 2.3   # right margin holds the pathogen labels
    ax_bar.set_ylabel("isolates", fontsize=7)
    ax_bar.set_xlim(-0.6, XR)
    ax_bar.set_ylim(0, max(c[1] for c in combos) * 1.18)
    ax_bar.set_xticks([])
    ax_bar.tick_params(labelsize=6, length=2)
    ax_bar.set_title("Pathogen specificity (taken forward)", fontsize=8, pad=4)

    # ---- membership matrix (pathogen labels on the RIGHT to clear panel B) ----
    for yi, p in enumerate(P):
        ax_mat.text(len(combos) - 0.35, yi, fs.patho_it(p), ha="left", va="center",
                    fontsize=6.5, color=fs.PATHO_COL[p])
        ax_mat.plot([-0.5, len(combos) - 0.5], [yi, yi], color="#eeeeee", lw=4, zorder=0)
    for xi, (cmb, _) in zip(x, combos):
        members = [yi for yi, p in enumerate(P) if p in cmb]
        for yi, p in enumerate(P):
            on = p in cmb
            ax_mat.scatter(xi, yi, s=34, color=fs.PATHO_COL[p] if on else "#d0d0d0",
                           edgecolor="none", zorder=3)
        if len(members) > 1:
            ax_mat.plot([xi, xi], [min(members), max(members)], color="#4a4a4a",
                        lw=1.4, zorder=2)
    ax_mat.set_xlim(-0.6, XR)
    ax_mat.set_ylim(len(P) - 0.5, -0.5)
    ax_mat.axis("off")


def species_legend(ax, d):
    from matplotlib.lines import Line2D
    present = [s for s in SP_ORDER if s in set(d.species)]
    handles = [Line2D([0], [0], marker="s", ls="", ms=6, color=spc(s), label=it(s))
               for s in present]
    if "Other" in set(d.species):
        handles.append(Line2D([0], [0], marker="s", ls="", ms=6, color="#bbbbbb",
                              label="sequencing failed"))
    ax.axis("off")
    ax.legend(handles=handles, loc="upper left", ncol=2, frameon=False,
              fontsize=6.5, handletextpad=0.3, columnspacing=1.0,
              labelspacing=0.6, borderaxespad=0)


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    d, strains = load_prelim()
    d = order_isolates(d, strains)
    m = pd.read_csv(f"{DATA}/phenotype_matrix.csv").set_index("isolate")
    fwd = sorted([i for i in d.isolate if i in m.index],
                 key=lambda i: -m.loc[i, "MI_Bcinerea"])

    fig = plt.figure(figsize=(7.2, 6.6))
    axA = fig.add_axes([0.055, 0.09, 0.28, 0.82])
    axB = fig.add_axes([0.46, 0.09, 0.16, 0.82])
    axC_bar = fig.add_axes([0.75, 0.62, 0.23, 0.28])
    axC_mat = fig.add_axes([0.75, 0.47, 0.23, 0.13])
    axLeg = fig.add_axes([0.72, 0.09, 0.28, 0.30])

    panel_A(axA, d, strains)
    panel_B(axB, fwd)
    panel_C(axC_bar, axC_mat)
    species_legend(axLeg, d)

    fs.panel_letter(axA, "A", x=-0.22, y=1.02)
    fs.panel_letter(axB, "B", x=-0.42, y=1.02)
    fs.panel_letter(axC_bar, "C", x=-0.18, y=1.06)

    fs.save(fig, f"{OUT}/Fig1_screen")
    print("wrote", f"{OUT}/Fig1_screen.svg/.png")
    print(f"panel A: {len(d)} isolates; panel B/C: {len(fwd)} taken forward")


if __name__ == "__main__":
    main()
