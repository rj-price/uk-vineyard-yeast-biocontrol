#!/usr/bin/env python
"""
20_fig2_phylo.py  ->  figures/manuscript/Fig2_phylogenomics.{svg,png}
Manuscript Figure 2 — phylogenomics & genome-based species identity.
(A) rectangular cladogram (BUSCO supermatrix, 55 isolates + type-strain references
    + Fusarium outgroup), tips coloured by genome-ANI species, clades bracketed.
(B) ANI heatmap for the taxonomically ambiguous genera (isolate x reference),
    yellow->red about the 95% species threshold.
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from Bio import Phylo
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs
from _palette import spc, it, SP_ORDER

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
DATA = f"{PROJ}/analysis/data"
SID = f"{PROJ}/analysis/species_id"
TREE = f"{PROJ}/analysis/phylo_refs/phylo_out/supermatrix/refs_tree.contree"
OUT = f"{PROJ}/analysis/figures/manuscript"
# sequenced but never screened (all Sporobolomyces, unphenotyped) -> excluded (56/55 study)
DROP = {"TY8", "TY9", "TY12", "TY17", "TY24"}


def norm_species(s):
    if not isinstance(s, str) or not s:
        return "Other"
    if s.startswith("Hanseniaspora"): return "Hanseniaspora uvarum"
    if s.startswith("Pichia"): return "Pichia sp."
    if s.startswith("Metschnikowia"): return "Metschnikowia sp."
    if "babjevae" in s: return "Rhodotorula babjevae"
    if s.startswith("Sporobolomyces roseus"): return "Sporobolomyces roseus"
    if s.startswith("Sporobolomyces"): return "Sporobolomyces sp."
    if "carnescens" in s: return "Vishniacozyma carnescens"
    if "victoriae" in s: return "Vishniacozyma victoriae"
    if s.startswith("Filobasidium"): return "Filobasidium sp."
    return "Other"


def sid(name):
    return re.sub(r"_EKDN.*", "", name)


# ---------- Panel A: cladogram ----------
def layout_cladogram(tree):
    """Right-aligned cladogram: x = -(edges to furthest descendant leaf),
    y = mean of children (leaves in ladderised order)."""
    leaves = tree.get_terminals()
    yy = {lf: i for i, lf in enumerate(leaves)}
    xx = {}

    def level(cl):
        if cl.is_terminal():
            xx[cl] = 0
            return 0
        h = 1 + max(level(c) for c in cl.clades)
        xx[cl] = -h
        return h
    level(tree.root)

    def ypos(cl):
        if cl.is_terminal():
            return yy[cl]
        y = np.mean([ypos(c) for c in cl.clades])
        yy[cl] = y
        return y
    ypos(tree.root)
    return xx, yy, leaves


def ref_label(name):
    """REF_Genus_species -> italic 'Genus species' with a superscript T (type strain)."""
    sp = name.replace("REF_", "").replace("_", " ")
    esc = sp.replace(" ", r"\ ")
    return f"$\\it{{{esc}}}^{{\\mathrm{{T}}}}$"


def panel_A(ax, iso2sp):
    tree = Phylo.read(TREE, "newick")
    # prune the sequenced-but-unscreened Sporobolomyces (56/55 study)
    for lf in list(tree.get_terminals()):
        if lf.name in DROP:
            tree.prune(lf)
    # root on the Fusarium outgroup
    outs = [l for l in tree.get_terminals() if l.name.startswith("OUT_")]
    if len(outs) >= 2:
        tree.root_with_outgroup(*outs)
    elif outs:
        tree.root_with_outgroup(outs[0])
    tree.ladderize()
    xx, yy, leaves = layout_cladogram(tree)

    def kind(name):
        if name.startswith("OUT_"):
            return "out"
        if name.startswith("REF_"):
            return "ref"
        return "iso"

    for cl in tree.get_nonterminals():
        cys = [yy[c] for c in cl.clades]
        ax.plot([xx[cl], xx[cl]], [min(cys), max(cys)], color="#555555", lw=0.6, zorder=1)
        for c in cl.clades:
            ax.plot([xx[cl], xx[c]], [yy[c], yy[c]], color="#555555", lw=0.6, zorder=1)
        if cl.confidence is not None and cl is not tree.root:
            s = float(cl.confidence)
            if s >= 95:
                ax.scatter(xx[cl], yy[cl], s=6, color="#333333", zorder=3, edgecolor="none")
            elif s >= 80:
                ax.scatter(xx[cl], yy[cl], s=6, facecolor="white",
                           edgecolor="#333333", lw=0.4, zorder=3)

    for lf in leaves:
        k = kind(lf.name)
        if k == "out":
            ax.scatter(0.15, yy[lf], s=11, color="#000000", edgecolor="none", zorder=4,
                       clip_on=False)
            ax.text(0.4, yy[lf], "$\\it{F.\\ venenatum}$", va="center", ha="left",
                    fontsize=4, color="#000000", clip_on=False)
        elif k == "ref":
            ax.scatter(0.15, yy[lf], s=13, marker="D", facecolor="white",
                       edgecolor="#888888", lw=0.6, zorder=4, clip_on=False)
            ax.text(0.4, yy[lf], ref_label(lf.name), va="center", ha="left", fontsize=4,
                    color="#888888", clip_on=False)
        else:
            sp = norm_species(iso2sp.get(lf.name, ""))
            ax.scatter(0.15, yy[lf], s=13, color=spc(sp), edgecolor="none", zorder=4,
                       clip_on=False)
            ax.text(0.4, yy[lf], lf.name, va="center", ha="left", fontsize=4,
                    color=fs.INK, clip_on=False)

    # species clade brackets (isolate tips only), black labels
    xb = 5.2
    groups = {}
    for lf in leaves:
        if kind(lf.name) == "iso":
            sp = norm_species(iso2sp.get(lf.name, ""))
            groups.setdefault(sp, []).append(yy[lf])
    for sp, ys in groups.items():
        if sp in ("Other",):
            continue
        y0, y1 = min(ys), max(ys)
        ax.plot([xb, xb], [y0, y1], color=spc(sp), lw=2.2, solid_capstyle="round",
                clip_on=False)
        ax.text(xb + 0.3, (y0 + y1) / 2, it(sp), va="center", ha="left", fontsize=5.5,
                color=fs.INK, clip_on=False)

    n_iso = sum(1 for lf in leaves if kind(lf.name) == "iso")
    ax.set_xlim(min(xx.values()) - 0.5, xb + 9)
    ax.set_ylim(len(leaves) - 0.5, -1.5)
    ax.axis("off")
    ax.scatter([], [], s=6, color="#333333", label="UFBoot ≥ 95")
    ax.scatter([], [], s=6, facecolor="white", edgecolor="#333333", lw=0.4,
               label="UFBoot 80–95")
    ax.scatter([], [], s=13, marker="D", facecolor="white", edgecolor="#888888", lw=0.6,
               label="reference genome$^{T}$")
    ax.legend(loc="lower left", frameon=False, fontsize=5.5, handletextpad=0.3,
              bbox_to_anchor=(0.0, 0.0))
    ax.set_title(f"BUSCO supermatrix phylogeny ({n_iso} isolates + references)",
                 fontsize=8, pad=2)


# ---------- Panel B: ANI heatmap ----------
GENERA = {
    "Rhodotorula cf. babjevae": ["Rhodotorula_babjevae", "Rhodotorula_graminis",
                                 "Rhodotorula_kratochvilovae", "Rhodotorula_mucilaginosa",
                                 "Rhodotorula_toruloides"],
    "Sporobolomyces": ["Sporobolomyces_roseus", "Rhodosporidiobolus_nylandii"],
    "Metschnikowia cf. pulcherrima": ["Metschnikowia_pulcherrima",
                                      "Metschnikowia_ziziphicola_fructicola",
                                      "Metschnikowia_chrysoperlae"],
    "Pichia kluyveri": ["Pichia_kluyveri", "Pichia_fermentans", "Pichia_membranifaciens"],
}


def panel_B(ax):
    a = pd.read_csv(f"{SID}/ani_raw.tsv", sep="\t", header=None,
                    names=["q", "r", "ani", "m", "t"])
    a["iso_sid"] = a.q.map(lambda p: sid(os.path.basename(p)))
    a["ref"] = a.r.map(lambda p: re.sub(r"__GC.*", "", os.path.basename(p)))
    m = pd.read_csv(f"{DATA}/isolate_master.csv")
    sid2iso = dict(zip(m.seq_id, m.isolate)); sid2sp = dict(zip(m.seq_id, m.species_genome))
    a["iso"] = a.iso_sid.map(sid2iso); a["sp"] = a.iso_sid.map(sid2sp)
    a = a[a.iso.notna()]

    # rows: isolates from the ambiguous genera, grouped; cols: their references
    genus_of = {"Rhodotorula cf. babjevae": lambda s: "babjevae" in str(s),
                "Sporobolomyces": lambda s: str(s).startswith("Sporobolomyces"),
                "Metschnikowia cf. pulcherrima": lambda s: str(s).startswith("Metschnikowia"),
                "Pichia kluyveri": lambda s: str(s).startswith("Pichia")}
    rows, row_lab, row_grpline = [], [], []
    cols = []
    for g, refs in GENERA.items():
        cols += refs
    seen = set(); cols = [c for c in cols if not (c in seen or seen.add(c))]

    order, group_of = [], []
    for g, test in genus_of.items():
        isos = sorted(i for i in a[a.sp.map(test)].iso.unique() if i not in DROP)
        for iso in isos:
            order.append(iso); group_of.append(g)
    M = np.full((len(order) + 4, len(cols)), np.nan)
    for ri, iso in enumerate(order):
        sub = a[a.iso == iso]
        for ci, c in enumerate(cols):
            v = sub[sub.ref == c].ani
            if len(v):
                M[ri, ci] = v.iloc[0]
    # Pichia is absent from ani_raw (fastANI reported only its best hit); pull the
    # kluyveri ANI from the per-isolate summary so all four ambiguous genera appear
    sa = pd.read_csv(f"{DATA}/species_id_ani.csv")
    pk = sa[sa.current_label.astype(str).str.startswith("Pichia")].sort_values("isolate")
    ci_k = cols.index("Pichia_kluyveri")
    base = len(order)
    for j, (_, r) in enumerate(pk.iterrows()):
        M[base + j, ci_k] = r.ani_best
        order.append(r.isolate); group_of.append("Pichia kluyveri")

    # drop reference columns with no data (e.g. Rhodosporidiobolus, distant Pichia refs)
    keep = ~np.all(np.isnan(M), axis=0)
    M = M[:, keep]
    cols = [c for c, k in zip(cols, keep) if k]

    from matplotlib.colors import Normalize
    norm = Normalize(vmin=78, vmax=100)   # yellow (low ANI) -> red (species-level)
    im = ax.imshow(M, aspect="auto", cmap="YlOrRd", norm=norm)
    im.cmap.set_bad("#f0f0f0")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("_", " ") for c in cols], rotation=40, ha="left",
                       fontsize=5.5, fontstyle="italic")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=4.8)
    for ri in range(len(order)):
        for ci in range(len(cols)):
            if not np.isnan(M[ri, ci]):
                ax.text(ci, ri, f"{M[ri,ci]:.0f}", ha="center", va="center", fontsize=4,
                        color="white" if M[ri, ci] >= 95 else fs.INK)
    # genus separators
    yb = 0
    for g in dict.fromkeys(group_of):
        n = group_of.count(g)
        if yb > 0:
            ax.axhline(yb - 0.5, color="white", lw=1.6)
        yb += n
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.03,
                      ticks=[80, 85, 90, 100], extend="min")
    cb.set_label("ANI (%)", fontsize=6)
    cb.ax.tick_params(labelsize=6, length=2)
    # mark the 95% species threshold on the colourbar
    cb.ax.axhline(95, color="#222222", lw=1.0)
    cb.ax.text(1.7, 95, "95 = species", fontsize=5, va="center", ha="left",
               color="#222222", transform=cb.ax.get_yaxis_transform())


# ---------- Panel C: corrections table ----------
def panel_C(ax):
    d = pd.read_csv(f"{DATA}/species_id_final.csv")
    ch = d[d.species_changed | d.genus_changed].copy()
    # collapse to distinct prior->new with counts
    g = (ch.groupby(["prior_label", "species_genome"])
         .agg(n=("isolate", "size"), genus=("genus_changed", "any")).reset_index()
         .sort_values("n", ascending=False))
    ax.axis("off")
    ax.set_title("ITS → genome corrections", fontsize=8, pad=4, loc="left")
    y = 0.98
    ax.text(0.02, y, "n", fontsize=6, fontweight="bold", transform=ax.transAxes)
    ax.text(0.10, y, "ITS / prior label", fontsize=6, fontweight="bold", transform=ax.transAxes)
    ax.text(0.55, y, "Genome ID", fontsize=6, fontweight="bold", transform=ax.transAxes)
    y -= 0.06
    ax.plot([0.02, 0.98], [y + 0.02, y + 0.02], color="#cccccc", lw=0.6,
            transform=ax.transAxes)
    for _, r in g.iterrows():
        spk = norm_species(r.species_genome)
        ax.text(0.02, y, str(r.n), fontsize=6, transform=ax.transAxes)
        tag = "  (genus)" if r.genus else ""
        ax.text(0.10, y, r.prior_label, fontsize=6, fontstyle="italic",
                transform=ax.transAxes, color=fs.MUTED)
        ax.text(0.55, y, r.species_genome + tag, fontsize=6, transform=ax.transAxes,
                color=spc(spk), fontweight="bold" if r.genus else "normal")
        y -= 0.075
    ax.text(0.02, y - 0.02, "Confirmed unchanged: Hanseniaspora uvarum, "
            "Vishniacozyma carnescens, Sporobolomyces roseus.",
            fontsize=5.2, transform=ax.transAxes, color=fs.MUTED, style="italic")


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    m = pd.read_csv(f"{DATA}/isolate_master.csv")
    iso2sp = dict(zip(m.isolate, m.species_genome))

    # ~78 tips (55 isolates + 21 refs + 2 outgroup) -> tall figure, two full-height columns
    fig = plt.figure(figsize=(7.4, 9.6))
    axA = fig.add_axes([0.01, 0.02, 0.45, 0.94])   # cladogram, full height, left
    axB = fig.add_axes([0.62, 0.30, 0.26, 0.56])   # ANI heatmap, top-right, labels on top

    panel_A(axA, iso2sp)
    panel_B(axB)
    fs.panel_letter(axA, "A", x=0.01, y=0.99)
    fs.panel_letter(axB, "B", x=-0.28, y=1.10)

    fs.save(fig, f"{OUT}/Fig2_phylogenomics")
    print("wrote", f"{OUT}/Fig2_phylogenomics.svg/.png")


if __name__ == "__main__":
    main()
