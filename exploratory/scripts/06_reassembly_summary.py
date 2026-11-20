#!/usr/bin/env python
"""
06_reassembly_summary.py
Compare re-assembled ascomycete genomes (analysis/reassembly/assemblies/<seq_id>/)
against the original SPAdes assemblies. Parses seqkit stats + BUSCO short summaries.
-> analysis/data/reassembly_comparison.csv  (+ prints a table)
"""
import os, re, glob
import pandas as pd, numpy as np

RA="analysis/reassembly"; DATA="analysis/data"
sheet=pd.read_csv(f"{RA}/samplesheet.tsv", sep="\t")

def busco_from_summary(path):
    if not path: return {}
    txt=open(path).read()
    m=re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", txt)
    if not m: return {}
    return dict(new_busco_C=float(m.group(1)), new_busco_S=float(m.group(2)),
                new_busco_D=float(m.group(3)), new_busco_F=float(m.group(4)),
                new_busco_M=float(m.group(5)))

def seqkit_stats(path):
    df=pd.read_csv(path, sep="\t")
    r=df.iloc[0]
    # seqkit -a columns include N50; find them case-insensitively
    cols={c.lower():c for c in df.columns}
    return dict(new_n_contigs=int(r[cols["num_seqs"]]), new_size=int(r[cols["sum_len"]]),
                new_n50=int(r[cols.get("n50","N50")]) if "n50" in cols or "N50" in df.columns else np.nan,
                new_max=int(r[cols["max_len"]]))

rows=[]
for _,s in sheet.iterrows():
    sid=s.seq_id; d=f"{RA}/assemblies/{sid}"
    st=f"{d}/{sid}.stats.tsv"
    if not os.path.exists(st): continue
    rec=dict(isolate=s.isolate, seq_id=sid, species=s.species, exp_mb=s.exp_genome_mb,
             old_size=s.old_size, old_n50=s.old_n50, old_busco=s.old_busco)
    rec.update(seqkit_stats(st))
    bs=glob.glob(f"{d}/short_summary.*.txt")
    if bs: rec.update(busco_from_summary(bs[0]))
    rows.append(rec)

if not rows:
    print("No completed re-assemblies found yet."); raise SystemExit
d=pd.DataFrame(rows)
d["size_ratio_old"]=(d.old_size/(d.exp_mb*1e6)).round(2)
d["size_ratio_new"]=(d.new_size/(d.exp_mb*1e6)).round(2)
d["busco_gain"]=(d.get("new_busco_C",np.nan)-d.old_busco).round(1)
d["n50_fold"]=(d.new_n50/d.old_n50).round(1)
cols=["isolate","seq_id","species","exp_mb","old_size","new_size","size_ratio_new",
      "old_n50","new_n50","n50_fold","old_busco","new_busco_C","new_busco_D","busco_gain","new_n_contigs"]
d=d[cols].sort_values(["species","seq_id"])
d.to_csv(f"{DATA}/reassembly_comparison.csv", index=False)
pd.set_option("display.width",220)
print(d.to_string(index=False))
print(f"\n{len(d)} re-assemblies summarised -> {DATA}/reassembly_comparison.csv")
print("\nMeans by species:")
print(d.groupby("species")[["new_size","new_n50","old_busco","new_busco_C"]].mean().round(0).to_string())
