#!/usr/bin/env python
"""
14_species_final.py
Consolidate genome-based species identification from whole-genome ANI (fastANI,
primary) + LSU D1/D2 barcode BLAST (for taxa with no conspecific reference genome).
Writes a final per-isolate table and adds a `species_genome` column to the master.
-> analysis/data/species_id_final.csv
"""
import os, re
import pandas as pd, numpy as np
SID="analysis/species_id"; DATA="analysis/data"
sub=pd.read_csv("analysis/mining/killer_subjects.tsv", sep="\t")
iso2cur=dict(zip(sub.isolate, sub.species.str.replace("_"," ")))
p2i={p:i for p,i in zip(sub.assembly,sub.isolate)}
for p,i in list(p2i.items()): p2i[os.path.abspath(p)]=i

ani=pd.read_csv(f"{SID}/ani_raw.tsv", sep="\t", header=None, names=["q","r","ani","m","t"])
ani["iso"]=ani.q.map(p2i)
ani["ref"]=ani.r.map(lambda p: re.sub(r"__GC.*","",os.path.basename(p)).replace("_"," "))
ani_best={iso:(g.sort_values('ani',ascending=False).iloc[0].ref,
               round(g.ani.max(),1)) for iso,g in ani.groupby("iso")}

# genome-based species calls, grouped by evidence (isolates -> final call)
G=lambda *x:list(x)
GROUPS=[
 # (isolates, final species, basis, confidence)
 (G(*[i for i in sub.isolate if iso2cur.get(i,"").startswith("Hanseniaspora")]),
   "Hanseniaspora uvarum","ANI 99.2%","species"),
 (G("RY12","RY13","RY18","RY23"),
   "Pichia kluyveri","ANI 99.0% (was 'P. fermentans' by ITS)","species"),
 (G("TY26","TY30"),
   "Metschnikowia pulcherrima","ANI 96.1% to M. pulcherrima type NRRL Y-7111 (near threshold; fuzzy species boundaries in clade); TY26&TY30 conspecific","species (provisional)"),
 (G("UY1","UY2","UY3","UY4","TY13","TY27","TY28","TY29"),
   "Rhodotorula aff. babjevae","ANI 90.3% to babjevae (>graminis 88.8%); sub-95% = distinct lineage; NOT R. graminis","near-species"),
 (G("TY1","TY2","TY3","TY6","TY15","TY16","TY18","TY20","TY22"),
   "Sporobolomyces roseus","ANI 96%","species"),
 (G("TY19","TY21","TY23","TY25"),
   "Sporobolomyces sp.","ANI 83% to S. roseus = distinct species (barcode: Sporidiobolales)","genus"),
 (G("TY4","TY5","TY11"),
   "Vishniacozyma carnescens","LSU D1/D2 100%","species (barcode)"),
 (G("TY7"),
   "Vishniacozyma victoriae","LSU D1/D2 100% (was 'Hannaella sp.' - GENUS CORRECTION)","species (barcode)"),
 (G("TY10"),
   "Vishniacozyma cf. victoriae","LSU D1/D2 99.8% (was 'Naganishia vaughanmartiniae' - GENUS CORRECTION)","species (barcode)"),
 (G("TY14"),
   "Filobasidium sp. (cf. globisporum)","LSU D1/D2 100% globisporum (was 'F. stepposum')","genus (barcode)"),
]
rows=[]
for isos,sp,basis,conf in GROUPS:
    for iso in isos:
        ab=ani_best.get(iso,("",np.nan))
        rows.append(dict(isolate=iso, prior_label=iso2cur.get(iso,""),
                         species_genome=sp, ani_best_ref=ab[0], ani=ab[1],
                         basis=basis, confidence=conf))
res=pd.DataFrame(rows)
res["genus_changed"]=res.prior_label.str.split().str[0]!=res.species_genome.str.split().str[0]
res["species_changed"]=(res.prior_label.str.replace(" sp.","",regex=False).str.strip()
                        !=res.species_genome.str.replace(r" (cf\.|sp\.).*","",regex=True).str.strip())
order={i:n for n,i in enumerate(sub.isolate)}
res=res.sort_values("isolate",key=lambda s:s.map(order)).reset_index(drop=True)
res.to_csv(f"{DATA}/species_id_final.csv", index=False)

m=pd.read_csv(f"{DATA}/isolate_master.csv")
m=m.drop(columns=[c for c in ("species_genome","id_basis") if c in m.columns])
m=m.merge(res[["isolate","species_genome","basis"]].rename(columns={"basis":"id_basis"}),
          on="isolate", how="left")
m.to_csv(f"{DATA}/isolate_master.csv", index=False)

pd.set_option("display.width",240); pd.set_option("display.max_rows",70)
print("=== final genome-based species (per group) ===")
print(res.drop_duplicates("species_genome")[["species_genome","basis","confidence"]].to_string(index=False))
print("\n=== isolates whose GENUS changed vs prior label ===")
print(res[res.genus_changed][["isolate","prior_label","species_genome"]].to_string(index=False))
print("\n=== species-level changes (same genus) ===")
sc=res[res.species_changed & ~res.genus_changed]
print(sc[["isolate","prior_label","species_genome"]].drop_duplicates("prior_label").to_string(index=False))
print(f"\nWrote {DATA}/species_id_final.csv + species_genome/id_basis in master ({len(res)} isolates)")
