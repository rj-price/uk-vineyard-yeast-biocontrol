#!/usr/bin/env python
"""
07_platanus_comparison.py
Compare original SPAdes ascomycete assemblies vs Platanus-allee re-assemblies
(heterozygosity-aware). Parses seqkit -a stats + BUSCO short summaries.
-> analysis/data/platanus_comparison.csv
   analysis/figures/fig10_platanus_vs_spades.png
"""
import os, re, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import SP_COL as sp_col, spc, it

RA="analysis/reassembly"; DATA="analysis/data"; FIG="analysis/figures"
master=pd.read_csv(f"{DATA}/isolate_master.csv")
sheet=pd.read_csv(f"{RA}/samplesheet.tsv", sep="\t")

def busco_pct(path):
    if not path: return {}
    m=re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", open(path).read())
    if not m: return {}
    return dict(new_busco=float(m.group(1)), new_S=float(m.group(2)), new_D=float(m.group(3)),
                new_F=float(m.group(4)), new_M=float(m.group(5)))

def sk_stats(path):
    df=pd.read_csv(path, sep="\t"); r=df.iloc[0]; c={x.lower():x for x in df.columns}
    return dict(new_contigs=int(r[c["num_seqs"]]), new_size=int(r[c["sum_len"]]),
                new_n50=int(r[c["n50"]]) if "n50" in c else np.nan,
                new_max=int(r[c["max_len"]]))

rows=[]
for _,s in sheet.iterrows():
    sid=s.seq_id; d=f"{RA}/platanus/{sid}"; st=f"{d}/{sid}.stats.tsv"
    if not os.path.exists(st): continue
    rec=dict(isolate=s.isolate, seq_id=sid, species=s.species, exp_mb=s.exp_genome_mb,
             old_size=s.old_size, old_n50=s.old_n50, old_busco=s.old_busco)
    rec.update(sk_stats(st))
    bs=glob.glob(f"{d}/short_summary.*.txt")
    if bs: rec.update(busco_pct(bs[0]))
    rows.append(rec)

if not rows:
    print("No Platanus assemblies found yet."); raise SystemExit
d=pd.DataFrame(rows)
d["size_ratio_new"]=(d.new_size/(d.exp_mb*1e6)).round(2)
d["busco_gain"]=(d.new_busco-d.old_busco).round(1)
d["n50_fold"]=(d.new_n50/d.old_n50).round(1)
cols=["isolate","seq_id","species","exp_mb","old_size","new_size","size_ratio_new",
      "old_n50","new_n50","n50_fold","old_busco","new_busco","new_D","busco_gain","new_contigs"]
d=d[cols].sort_values(["species","new_busco"], ascending=[True,False])
d.to_csv(f"{DATA}/platanus_comparison.csv", index=False)
pd.set_option("display.width",240)
print(d.to_string(index=False))
print(f"\n{len(d)} isolates. Mean BUSCO: SPAdes {d.old_busco.mean():.1f}% -> Platanus {d.new_busco.mean():.1f}% "
      f"(+{d.busco_gain.mean():.1f}); mean N50 {d.old_n50.mean():.0f} -> {d.new_n50.mean():.0f} bp")

# ---- FIGURE 10: paired BUSCO before/after, + N50 ----
fig, axes = plt.subplots(1,2, figsize=(14,6))
d2=d.sort_values(["species","new_busco"]).reset_index(drop=True)
y=np.arange(len(d2))
ax=axes[0]
for i,r in d2.iterrows():
    ax.plot([r.old_busco,r.new_busco],[i,i],color="#ccc",lw=1.4,zorder=1)
ax.scatter(d2.old_busco,y,s=45,color="#bbbbbb",edgecolor="white",zorder=2,label="SPAdes (original)")
ax.scatter(d2.new_busco,y,s=55,color=[spc(s) for s in d2.species],edgecolor="white",zorder=3,label="Platanus-allee")
ax.set_yticks(y); ax.set_yticklabels(d2.isolate,fontsize=7)
ax.set_xlabel("BUSCO completeness (%, fungi_odb10)"); ax.set_xlim(0,100)
ax.axvline(90,color="#009E73",ls="--",lw=0.8); ax.set_title("BUSCO completeness: original → re-assembly",fontsize=11,loc="left")
sp_present=[s for s in sp_col if s in set(d2.species)]
h=[Patch(facecolor="#bbbbbb",label="SPAdes (original)")]+[Patch(facecolor=sp_col[s],label=it(s)) for s in sp_present]
ax.legend(handles=h,frameon=False,fontsize=8,loc="lower right")

ax=axes[1]
ax.barh(y-0.2,d2.old_n50/1000,height=0.4,color="#bbbbbb",label="SPAdes")
ax.barh(y+0.2,d2.new_n50/1000,height=0.4,color=[spc(s) for s in d2.species],label="Platanus")
ax.set_yticks(y); ax.set_yticklabels(d2.isolate,fontsize=7)
ax.set_xlabel("Scaffold N50 (kb)"); ax.set_title("Contiguity (N50): original → re-assembly",fontsize=11,loc="left")
fig.suptitle("Platanus-allee rescues the heterozygous ascomycete genomes that SPAdes could not assemble",
             x=0.5,y=1.01,fontsize=13)
plt.savefig(f"{FIG}/fig10_platanus_vs_spades.png",dpi=300,bbox_inches="tight"); plt.close()
print(f"\nWrote {DATA}/platanus_comparison.csv + {FIG}/fig10_platanus_vs_spades.png")
