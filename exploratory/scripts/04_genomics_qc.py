#!/usr/bin/env python
"""
04_genomics_qc.py
Assembly QC summary + figure, and kraken2 classification-rate summary.
-> analysis/figures/fig7_assembly_qc.png
   analysis/data/assembly_qc_summary.csv
   analysis/data/kraken2_summary.csv
"""
import os, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import SP_ORDER as sp_all, SP_COL as sp_col, it

DATA="analysis/data"; FIG="analysis/figures"
master = pd.read_csv(f"{DATA}/isolate_master.csv")

# ---- QC quality tier ----
d = master.dropna(subset=["busco_complete","n50"]).copy()
def tier(r):
    if r.busco_complete>=90 and r.n50>=100000: return "High"
    if r.busco_complete>=80 and r.n50>=50000:  return "Medium"
    return "Low (fragmented)"
d["qc_tier"] = d.apply(tier, axis=1)

qc = d[["isolate","genome_code","seq_id","species","size_inclN","gc","n50","seq_num",
        "busco_complete","busco_dup","busco_frag","busco_missing","qc_tier"]].copy()
qc["assembly_Mb"] = (qc["size_inclN"]/1e6).round(2)
qc = qc.sort_values(["qc_tier","species","isolate"])
qc.to_csv(f"{DATA}/assembly_qc_summary.csv", index=False)

# per-species QC means
sp_qc = d.groupby("species").agg(
    n=("isolate","size"), assembly_Mb=("size_inclN", lambda x:(x.mean()/1e6)),
    GC=("gc","mean"), N50_kb=("n50", lambda x:x.mean()/1e3),
    BUSCO_complete=("busco_complete","mean"), BUSCO_missing=("busco_missing","mean")
).round(2).sort_values("BUSCO_complete", ascending=False)
sp_qc.to_csv(f"{DATA}/assembly_qc_by_species.csv")

# ---- FIGURE 7: BUSCO completeness vs N50, coloured by species, sized by assembly ----
fig, ax = plt.subplots(figsize=(8,5.6))
for s in sp_all:
    sub = d[d.species==s]
    if len(sub):
        ax.scatter(sub.n50/1e3, sub.busco_complete, s=(sub.size_inclN/1e6)*7,
                   color=sp_col[s], edgecolor="white", linewidth=0.7, alpha=0.9,
                   label=it(s), zorder=3)
ax.set_xscale("log")
ax.axhspan(90,100, color="#009E73", alpha=0.06)
ax.axhline(90, color="#009E73", lw=0.8, ls="--")
ax.axhline(80, color="#E69F00", lw=0.8, ls="--")
ax.text(ax.get_xlim()[1], 90, " 90% complete", color="#009E73", fontsize=7, va="bottom", ha="right")
ax.set_xlabel("Scaffold N50 (kb, log scale)")
ax.set_ylabel("BUSCO completeness (%, fungi_odb10)")
ax.set_title("Assembly quality splits basidiomycete from ascomycete yeasts",
             fontsize=11, loc="left")
ax.legend(frameon=False, fontsize=7.5, title="Species", title_fontsize=8,
          loc="lower right", ncol=1)
ax.text(0.02,0.03,"marker size ∝ assembly size", transform=ax.transAxes,
        fontsize=7, color="#777", style="italic")
plt.savefig(f"{FIG}/fig7_assembly_qc.png", dpi=300, bbox_inches="tight"); plt.close()

# ---- kraken2 summary (job-level; logs carry no sample id) ----
rows=[]
for f in sorted(glob.glob("genomics/kraken2/slurm-*.out")):
    txt=open(f).read()
    m=re.search(r"([\d.]+)%\)\s*\n\s*\d+ sequences unclassified", txt) or \
      re.search(r"(\d+) sequences classified \(([\d.]+)%\)", txt)
    seqs=re.search(r"(\d+) sequences \(([\d.]+) Mbp\)", txt)
    cl=re.search(r"(\d+) sequences classified \(([\d.]+)%\)", txt)
    if cl:
        rows.append((os.path.basename(f), int(seqs.group(1)) if seqs else np.nan,
                     float(seqs.group(2)) if seqs else np.nan,
                     int(cl.group(1)), float(cl.group(2))))
k=pd.DataFrame(rows, columns=["logfile","n_seqs","Mbp","n_classified","pct_classified"])
k.to_csv(f"{DATA}/kraken2_summary.csv", index=False)

print("=== Assembly QC tiers ===")
print(qc.qc_tier.value_counts().to_string())
print("\n=== QC by species ===")
print(sp_qc.to_string())
print(f"\n=== kraken2 ({len(k)} assemblies) classified % ===")
print(f"median {k.pct_classified.median():.1f}%  range {k.pct_classified.min():.1f}-{k.pct_classified.max():.1f}%")
print("Wrote fig7 + QC/kraken2 summary tables")
