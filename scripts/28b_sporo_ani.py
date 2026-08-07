#!/usr/bin/env python
"""
28b_sporo_ani.py -> Supplementary figure: ANI of the Sporobolomyces isolates against
ALL Sporobolomyces reference species (from the expanded fastANI run ani2_raw.tsv).
Shows the S. roseus group (>=95% to S. roseus) versus the cryptic second species
(~83% to S. roseus, <=78% to every other Sporobolomyces species = undescribed).
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SID = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id"
FIG = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/figures/manuscript"
os.makedirs(FIG, exist_ok=True)

df = pd.read_csv(f"{SID}/ani2_raw.tsv", sep="\t", header=None,
                 names=["q", "r", "ani", "frag", "total"])
df["q"] = df["q"].map(lambda p: os.path.basename(p).replace(".fasta", ""))
df["rsp"] = df["r"].map(lambda p: os.path.basename(p).replace(".fna", "").split("__")[0])

# Sporobolomyces isolates = those whose best hit is a Sporobolomyces reference
best = df.sort_values("ani", ascending=False).groupby("q").head(1)
iso_sporo = sorted(best[best.rsp.str.startswith("Sporobolomyces")].q,
                   key=lambda s: (s[:2], int("".join(ch for ch in s if ch.isdigit()) or 0)))

# full set of Sporobolomyces references on disk (some are too divergent for fastANI
# to report a value for these isolates; those columns stay blank = below detection)
import glob
sporo_refs = sorted({os.path.basename(f).split("__")[0]
                     for f in glob.glob(f"{SID}/refs/Sporobolomyces_*.fna")})
# column order: roseus, pararoseus, salmonicolor (the three kept in the main figure) first
lead = ["Sporobolomyces_roseus", "Sporobolomyces_pararoseus", "Sporobolomyces_salmonicolor"]
cols = [c for c in lead if c in sporo_refs] + [c for c in sporo_refs if c not in lead]

M = (df[df.q.isin(iso_sporo) & df.rsp.isin(cols)]
     .pivot_table(index="q", columns="rsp", values="ani", aggfunc="max")
     .reindex(index=iso_sporo, columns=cols))

# split isolates into roseus-group (>=95% to S. roseus) and cryptic (<95%)
ros = M["Sporobolomyces_roseus"]
roseus_grp = [i for i in iso_sporo if ros.get(i, 0) >= 95]
cryptic = [i for i in iso_sporo if ros.get(i, 0) < 95]
row_order = roseus_grp + cryptic
M = M.reindex(index=row_order)

def sp_lab(r):
    g, s = r.split("_")[0], r.split("_")[1]
    return f"$\\it{{{g[0]}.\\ {s}}}$"

fig, ax = plt.subplots(figsize=(7.4, 5.4))
data = M.values.astype(float)
im = ax.imshow(data, cmap="YlGnBu", vmin=76, vmax=100, aspect="auto")

ax.set_xticks(range(len(cols))); ax.set_xticklabels([sp_lab(c) for c in cols], rotation=40, ha="right", fontsize=8)
ax.set_yticks(range(len(row_order))); ax.set_yticklabels(row_order, fontsize=8)

# annotate cells
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        v = data[i, j]
        if np.isnan(v):
            ax.text(j, i, "–", ha="center", va="center", fontsize=7, color="#999999")
        else:
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7,
                    color="white" if v >= 92 else "black")
            if v >= 95:  # species-level match: outline
                ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=False, ec="#d62728", lw=2))

# separator between roseus group and cryptic species
sep = len(roseus_grp) - 0.5
ax.axhline(sep, color="black", lw=1.4)
ax.annotate("$\\it{S.\\ roseus}$ group", xy=(-0.62, (len(roseus_grp)-1)/2),
            xycoords=("axes fraction", "data"), rotation=90, va="center", ha="center", fontsize=8)
ax.annotate("cryptic sp.", xy=(-0.62, (len(roseus_grp)+len(cryptic)-1+len(roseus_grp))/2 - 0.0),
            xycoords=("axes fraction", "data"), rotation=90, va="center", ha="center", fontsize=8)

cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
cb.set_label("ANI (%)", fontsize=8); cb.ax.tick_params(labelsize=7)
ax.set_title("Whole-genome ANI of the $\\it{Sporobolomyces}$ isolates against all\n"
             "sampled $\\it{Sporobolomyces}$ species (red outline = ≥95% species threshold)",
             fontsize=9)
fig.text(0.5, -0.02,
         "Blank cells (–): ANI below fastANI's reporting threshold (~76%), i.e. more "
         "divergent than the values shown.",
         ha="center", va="top", fontsize=7, color="#444444")
plt.tight_layout()
for ext in ("png", "svg"):
    plt.savefig(f"{FIG}/FigS4_sporobolomyces_ani.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("wrote FigS4_sporobolomyces_ani.png/.svg")
print(f"isolates: {len(roseus_grp)} roseus-group + {len(cryptic)} cryptic; refs: {len(cols)}")
print("cryptic max-ANI to any Sporobolomyces ref:",
      {i: round(np.nanmax(M.loc[i].values), 1) for i in cryptic})
