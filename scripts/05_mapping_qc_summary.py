#!/usr/bin/env python
"""
05_mapping_qc_summary.py
Parse fastp JSON + samtools flagstat/coverage for each mapped sample into one
table, compute coverage vs expected genome size and vs the assembly, and plot.
-> analysis/data/mapping_qc.csv ; analysis/figures/fig8_mapping_qc.png
"""
import os, json, re, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import SP_COL as sp_col, spc, it

QC="analysis/mapping_qc"; RES=f"{QC}/results"; DATA="analysis/data"; FIG="analysis/figures"
sheet=pd.read_csv(f"{QC}/samplesheet.tsv", sep="\t")

def parse_flagstat(p):
    d={}
    for line in open(p):
        if "in total" in line: d["total"]=int(line.split()[0])
        elif re.search(r"\bmapped \(", line): d["mapped"]=int(line.split()[0]); \
             d["pct_mapped"]=float(re.search(r"\(([\d.]+)%",line).group(1))
        elif "properly paired" in line:
            d["pct_proper"]=float(re.search(r"\(([\d.]+)%",line).group(1))
        elif "with mate mapped to a different chr" in line and "mapQ" not in line:
            d["mate_diff_chr"]=int(line.split()[0])
        elif "duplicates" in line: d["dups"]=int(line.split()[0])
    return d

def parse_coverage(p):  # samtools coverage: weighted mean depth + breadth over contigs
    df=pd.read_csv(p, sep="\t")
    L=(df["endpos"]-df["startpos"]+1)
    tot=L.sum()
    return dict(mean_depth=float((df["meandepth"]*L).sum()/tot),
                breadth_pct=float(100*df["covbases"].sum()/tot),
                n_contigs=len(df))

rows=[]
for _,r in sheet.iterrows():
    sid=r.seq_id
    fj=f"{RES}/{sid}.fastp.json"; ff=f"{RES}/{sid}.flagstat.txt"; fc=f"{RES}/{sid}.coverage.txt"
    if not (os.path.exists(fj) and os.path.exists(ff) and os.path.exists(fc)):
        print("MISSING outputs for", sid); continue
    j=json.load(open(fj))
    bb=j["summary"]["before_filtering"]; af=j["summary"]["after_filtering"]
    insert_peak=j.get("insert_size",{}).get("peak",np.nan)
    dup=j.get("duplication",{}).get("rate",np.nan)
    adapter=j.get("adapter_cutting",{}).get("adapter_trimmed_reads",0)
    fs=parse_flagstat(ff); cov=parse_coverage(fc)
    raw_bases=bb["total_bases"]
    rows.append(dict(
        isolate=r.isolate, seq_id=sid, species=r.species, group=r.group,
        reads_M=round(bb["total_reads"]/1e6,1),
        read_len=bb.get("read1_mean_length",np.nan),
        raw_Gb=round(raw_bases/1e9,2),
        q30_pct=round(100*bb["q30_rate"],1) if bb.get("q30_rate") else round(100*bb["q30_bases"]/raw_bases,1),
        gc_pct=round(100*bb["gc_content"],1),
        dup_pct=round(100*dup,1) if dup==dup else np.nan,
        adapter_pct=round(100*adapter/bb["total_reads"],1),
        insert_peak=insert_peak,
        cov_vs_genome=round(raw_bases/(r.exp_genome_mb*1e6),0),
        cov_vs_assembly=round(raw_bases/r.asm_size,0),
        pct_mapped=fs.get("pct_mapped",np.nan),
        pct_proper=fs.get("pct_proper",np.nan),
        pct_diff_contig=round(100*fs.get("mate_diff_chr",0)/fs["total"],1) if fs.get("total") else np.nan,
        mean_depth=round(cov["mean_depth"],0),
        breadth_pct=round(cov["breadth_pct"],1),
        busco=r.busco, n50=r.n50,
    ))
d=pd.DataFrame(rows).sort_values(["group","species","isolate"])
d.to_csv(f"{DATA}/mapping_qc.csv", index=False)
pd.set_option("display.width",200)
print(d.to_string(index=False))

# ---- FIGURE 8: mapping QC, ascomycetes vs basidiomycete controls ----
fig, axes = plt.subplots(1,4, figsize=(16,5))
d=d.reset_index(drop=True)
cols=[spc(s) for s in d.species]
def bars(ax, y, title, ylab, hline=None):
    ax.bar(range(len(d)), d[y], color=cols, edgecolor="white")
    ax.set_xticks(range(len(d))); ax.set_xticklabels(d.isolate, rotation=90, fontsize=7)
    ax.set_title(title, fontsize=10.5, loc="left"); ax.set_ylabel(ylab)
    nc=(d.group=="ascomycete").sum()
    ax.axvline(nc-0.5, color="#888", ls="--", lw=0.8)   # ascomycete | control divider
    if hline is not None: ax.axhline(hline, color="#bbb", ls=":", lw=0.8)
bars(axes[0],"pct_mapped","Reads map to own assembly","% reads mapped")
bars(axes[1],"pct_diff_contig","Mate on a DIFFERENT contig\n(fragmentation signature)","% of reads")
bars(axes[2],"cov_vs_genome","Raw coverage vs expected genome","fold coverage (x)")
bars(axes[3],"dup_pct","Duplication (fastp)","% duplicated reads")
from matplotlib.patches import Patch
present=[s for s in sp_col if s in set(d.species)]
handles=[Patch(facecolor=sp_col[s], label=it(s)) for s in present]
axes[0].legend(handles=handles, frameon=False, fontsize=7.5, loc="lower left")
fig.suptitle("Raw-read mapping QC: ascomycetes (left of dashed line) vs basidiomycete controls — "
             "reads are fine; the assemblies are shattered", x=0.5, y=1.02, fontsize=12)
plt.savefig(f"{FIG}/fig8_mapping_qc.png", dpi=300, bbox_inches="tight"); plt.close()
print("\nWrote mapping_qc.csv + fig8_mapping_qc.png")
