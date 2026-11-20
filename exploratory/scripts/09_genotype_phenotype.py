#!/usr/bin/env python
"""
09_genotype_phenotype.py
Cross-reference mined biocontrol genes against Gideon's phenotype assays for the
isolates that were both sequenced+mined AND assayed.
-> analysis/data/genotype_phenotype.csv  (+ printed narrative table)
"""
import os, re
import numpy as np, pandas as pd
RAW="analysis/data/raw_csv"; DATA="analysis/data"
inv=pd.read_csv(f"{DATA}/mining_gene_inventory.csv")
long=pd.read_csv(f"{DATA}/phenotypes_long.csv")

# qualitative lytic-enzyme assay (majority of 3 reps == positive)
q=pd.read_csv(f"{RAW}/Lytic_enzyme_assay_results__Qualitative.csv", header=None, dtype=str).fillna("")
enz_names=["Cellulase","Amylase","Lipase","Chitinase","Protease","Esterase","Beta-glucosidase"]
qhdr_row=2  # data starts row index 2 (0-based) after two header rows
qual={}
for _,r in q.iloc[2:].iterrows():
    iso=str(r[0]).strip()
    if not re.match(r"^(UY|TY|RY)\d+$", iso): continue
    calls={}
    for ei,en in enumerate(enz_names):
        c0=1+ei*3
        yes=sum(1 for k in range(3) if str(r[c0+k]).strip().lower()=="yes")
        calls[en]= "＋" if yes>=2 else "－"
    qual[iso]=calls

# %MI helpers
def mi(assay,iso,pathogen=None,cond=None):
    s=long[(long.assay==assay)&(long.isolate==iso)]
    if pathogen: s=s[s.pathogen==pathogen]
    if cond: s=s[s.condition==cond]
    return round(s.pct_MI.mean(),1) if len(s) else np.nan

rows=[]
for _,g in inv.iterrows():
    iso=g.isolate
    qc=qual.get(iso,{})
    rows.append(dict(
        isolate=iso, species=g.species,
        # --- lytic enzymes: genotype (gene count) vs phenotype (assay call) ---
        chitinase_genes=int(g["Chitinase (GH18/19)"]), chitinase_assay=qc.get("Chitinase","n/a"),
        glucanase_genes=int(g["b-1,3-glucanase (GH16/17/55/64/81/128)"]), bglucosidase_assay=qc.get("Beta-glucosidase","n/a"),
        cellulase_genes=int(g["Cellulase (GH5/6/7/12/45)"]), cellulase_assay=qc.get("Cellulase","n/a"),
        amylase_genes=int(g["Amylase (GH13/15)"]), amylase_assay=qc.get("Amylase","n/a"),
        protease_genes=int(g["MEROPS_protease"]), protease_assay=qc.get("Protease","n/a"),
        # --- iron competition: pulcherrimin/siderophore genotype vs iron %MI ---
        pulcherrimin_CDPS=int(g.get("SM:CDPS / pulcherrimin (iron)",0)),
        iron_MI_Bc=mi("iron",iso,"B. cinerea"),
        # --- VOCs: terpene clusters vs VOC %MI ---
        terpene_clusters=int(g.get("SM:Terpenes (VOCs)",0)),
        voc_MI_Bc=mi("VOC",iso,"B. cinerea"),
        # --- overall antagonism ---
        dual_MI_Bc=mi("dual",iso,"B. cinerea"),
        dual_MI_Fo=mi("dual",iso,"F. oxysporum"),
    ))
gp=pd.DataFrame(rows)
gp.to_csv(f"{DATA}/genotype_phenotype.csv", index=False)

pd.set_option("display.width",260); pd.set_option("display.max_columns",50)
print("=== GENOTYPE (gene counts) vs PHENOTYPE (assay) ===\n")
print("Lytic enzymes — gene count | assay call (＋/－):")
for _,r in gp.iterrows():
    print(f"  {r.isolate:5s} {r.species:26s} "
          f"chit {r.chitinase_genes:>2}|{r.chitinase_assay}  "
          f"gluc {r.glucanase_genes:>2}|{r.bglucosidase_assay}  "
          f"cell {r.cellulase_genes:>2}|{r.cellulase_assay}  "
          f"amyl {r.amylase_genes:>2}|{r.amylase_assay}  "
          f"prot {r.protease_genes:>3}|{r.protease_assay}")
print("\nIron competition & VOCs — cluster genotype vs %MI phenotype:")
for _,r in gp.iterrows():
    print(f"  {r.isolate:5s} {r.species:26s} "
          f"pulcherrimin_CDPS={r.pulcherrimin_CDPS}  iron%MI(Bc)={r.iron_MI_Bc}   "
          f"terpene_clusters={r.terpene_clusters}  VOC%MI(Bc)={r.voc_MI_Bc}")
print(f"\nWrote {DATA}/genotype_phenotype.csv")
