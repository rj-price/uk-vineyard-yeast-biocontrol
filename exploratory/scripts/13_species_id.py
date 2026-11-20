#!/usr/bin/env python
"""
13_species_id.py
Genome-based species identification for all 60 isolates via whole-genome ANI
(fastANI vs 19 type/representative references). Assigns species by best ANI:
  >=95%  confident species-level match
  90-95% same genus, species uncertain (or no conspecific reference available)
  <90%   genus-level only / no close reference
Compares against the current master labels and the ITS-based calls, flags changes.
-> analysis/data/species_id_ani.csv
"""
import os, re
import pandas as pd, numpy as np
SID="analysis/species_id"; DATA="analysis/data"
ani=pd.read_csv(f"{SID}/ani_raw.tsv", sep="\t", header=None,
                names=["q","r","ani","matched","total"])
sub=pd.read_csv("analysis/mining/killer_subjects.tsv", sep="\t")
path2iso=dict(zip(sub.assembly, sub.isolate))
iso2cur=dict(zip(sub.isolate, sub.species))
ani["isolate"]=ani.q.map(lambda p: path2iso.get(p) or path2iso.get(os.path.abspath(p)))
ani["ref_species"]=ani.r.map(lambda p: re.sub(r"__GC[AF].*","",os.path.basename(p)))
ani["aln_frac"]=(ani.matched/ani.total).round(3)

rows=[]
for iso,g in ani.groupby("isolate"):
    g=g.sort_values("ani",ascending=False)
    best=g.iloc[0]
    second=g.iloc[1] if len(g)>1 else None
    call=("species" if best.ani>=95 else "genus" if best.ani>=90 else "ambiguous")
    rows.append(dict(
        isolate=iso, current_label=iso2cur.get(iso,"").replace("_"," "),
        ani_best_species=best.ref_species.replace("_"," "), ani_best=round(best.ani,2),
        aln_frac=best.aln_frac, call=call,
        second_species=second.ref_species.replace("_"," ") if second is not None else "",
        second_ani=round(second.ani,2) if second is not None else np.nan,
    ))
res=pd.DataFrame(rows)
# normalise for comparison (genus-level agreement / species change)
def genus(s): return str(s).split()[0] if s else ""
res["genus_current"]=res.current_label.map(genus)
res["genus_ani"]=res.ani_best_species.map(genus)
res["changed"]=(res.current_label.str.replace(" sp.","",regex=False).str.strip()
                != res.ani_best_species.str.strip()) & (res.call=="species")
res=res.sort_values(["genus_ani","isolate"]).reset_index(drop=True)
res.to_csv(f"{DATA}/species_id_ani.csv", index=False)

pd.set_option("display.width",220); pd.set_option("display.max_rows",70)
print(res[["isolate","current_label","ani_best_species","ani_best","aln_frac","call",
           "second_species","second_ani"]].to_string(index=False))
print("\n--- species-level reassignments vs current label ---")
ch=res[res.changed & (res.genus_current==res.genus_ani)]
print(ch[["isolate","current_label","ani_best_species","ani_best"]].to_string(index=False) if len(ch) else "(none within same genus)")
print("\n--- genus disagreements (check!) ---")
gd=res[res.genus_current!=res.genus_ani]
print(gd[["isolate","current_label","ani_best_species","ani_best","call"]].to_string(index=False) if len(gd) else "(none)")
print("\n--- low-confidence (no conspecific reference; <95% best ANI) ---")
lc=res[res.call!="species"]
print(lc[["isolate","current_label","ani_best_species","ani_best","call"]].to_string(index=False) if len(lc) else "(none)")
print(f"\nWrote {DATA}/species_id_ani.csv")
