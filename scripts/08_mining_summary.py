#!/usr/bin/env python
"""
08_mining_summary.py
Summarise biocontrol-relevant genes across the mined genomes and cross-reference
with Gideon's phenotype assays.
Parses per genome (analysis/mining/results/<iso>/):
  dbcan/overview.tsv          -> CAZyme families (consensus >=2 tools) -> enzyme classes
  <iso>.annotations.txt       -> MEROPS proteases, antifungal/killer-toxin candidates
  antismash/<iso>.json        -> secondary-metabolite cluster types
Outputs:
  analysis/data/mining_gene_inventory.csv
  analysis/figures/fig11_biocontrol_genes.png
"""
import os, re, json, glob, collections
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import SP_COL as sp_col, spc, it

MINE="analysis/mining"; DATA="analysis/data"; FIG="analysis/figures"
sheet=pd.read_csv(f"{MINE}/samplesheet.tsv", sep="\t")

# CAZy family -> biocontrol enzyme class (maps to Gideon's lytic-enzyme assays)
CAZY_CLASS = {
 "Chitinase (GH18/19)":       ["GH18","GH19"],
 "Chitin-binding (CBM18/50)": ["CBM18","CBM50"],
 "b-1,3-glucanase (GH16/17/55/64/81/128)": ["GH16","GH17","GH55","GH64","GH81","GH128","GH152"],
 "b-1,6-glucanase (GH30)":    ["GH30"],
 "Cellulase (GH5/6/7/12/45)": ["GH5","GH6","GH7","GH12","GH45"],
 "Amylase (GH13/15)":         ["GH13","GH15"],
 "Mannanase (GH76)":          ["GH76"],
 "Pectinase (GH28/PL1/PL3/PL4)": ["GH28","PL1","PL3","PL4"],
}
# antiSMASH product -> mechanism bucket
SM_BUCKET = {
 "terpene":"Terpenes (VOCs)","terpene-precursor":"Terpenes (VOCs)",
 "fungal_CDPS":"CDPS / pulcherrimin (iron)","NI-siderophore":"Siderophore (iron)",
 "siderophore":"Siderophore (iron)",
 "NRPS":"NRPS","NRPS-like":"NRPS","NRPS-independent-siderophore":"Siderophore (iron)",
 "T1PKS":"PKS","T3PKS":"PKS","PKS-like":"PKS","betalactone":"Other","indole":"Other",
}

def dbcan_families(path):
    df=pd.read_csv(path, sep="\t"); df.columns=[c.strip() for c in df.columns]
    df["n"]=pd.to_numeric(df["#ofTools"],errors="coerce").fillna(0)
    cons=df[df["n"]>=2]
    fam=collections.Counter()
    for _,r in cons.iterrows():
        for col in ("dbCAN_hmm","DIAMOND","dbCAN_sub"):
            m=re.search(r"(GH\d+|GT\d+|CBM\d+|CE\d+|PL\d+|AA\d+)", str(r[col]))
            if m: fam[m.group(1)]+=1; break
    return fam, len(cons)

def antismash_products(iso, d):
    js=glob.glob(f"{d}/antismash/*.json")
    prods=collections.Counter()
    if js:
        data=json.load(open(js[0]))
        for rec in data.get("records",[]):
            for a in rec.get("areas",[]):
                for p in a.get("products",[]): prods[p]+=1
    return prods

rows=[]; sm_all=collections.Counter()
for _,s in sheet.iterrows():
    iso=s.isolate; d=f"{MINE}/results/{iso}"
    ov=f"{d}/dbcan/overview.tsv"; ann=f"{d}/{iso}.annotations.txt"
    if not os.path.exists(ov):
        print("skip (no dbcan yet):", iso); continue
    fam,ncaz=dbcan_families(ov)
    rec=dict(isolate=iso, species=s.species.replace("_"," "), n_cazyme=ncaz)
    for cls,fl in CAZY_CLASS.items():
        rec[cls]=sum(fam.get(f,0) for f in fl)
    # proteases + antifungal candidates from annotations
    if os.path.exists(ann):
        a=pd.read_csv(ann,sep="\t",dtype=str).fillna("")
        rec["MEROPS_protease"]=(a["Protease"]!="").sum()
        rec["n_genes"]=len(a)
        rec["antifungal_named"]=a["Product"].str.contains(
            "killer|antifungal|zymocin|pulcherrimin",case=False,na=False).sum()
    # antiSMASH clusters -> mechanism buckets
    prods=antismash_products(iso,d)
    smb=collections.Counter()
    for p,c in prods.items(): smb[SM_BUCKET.get(p,"Other")]+=c
    for b,c in smb.items(): rec[f"SM:{b}"]=c
    rec["n_SM_clusters"]=sum(prods.values())
    sm_all.update(smb)
    rows.append(rec)

if not rows:
    print("No mined genomes ready yet."); raise SystemExit
inv=pd.DataFrame(rows).fillna(0)
# species order for display
inv=inv.sort_values("species").reset_index(drop=True)
inv.to_csv(f"{DATA}/mining_gene_inventory.csv", index=False)
pd.set_option("display.width",240)
show=["isolate","species","n_genes","n_cazyme","MEROPS_protease","n_SM_clusters"]+list(CAZY_CLASS)
print(inv[[c for c in show if c in inv]].to_string(index=False))
print("\nSecondary-metabolite clusters across genomes:", dict(sm_all))

# ---- FIGURE 11: biocontrol gene inventory heatmap (genomes x categories) ----
# order genomes: Ascomycota (Hanseniaspora, Pichia, Metschnikowia) then Basidiomycota
order=["RY14","RY13","TY26","TY30","UY1","TY6","TY4"]
inv=inv.set_index("isolate").loc[[i for i in order if i in inv.isolate.values]].reset_index()
cat_cols=list(CAZY_CLASS)+["MEROPS_protease"]+[c for c in inv.columns if c.startswith("SM:")]
H=inv.set_index("isolate")[cat_cols].astype(float)
H=H.loc[:, (H.sum(axis=0)>0)]              # drop all-zero categories
cat_cols=list(H.columns)
labels=[f"{i}\n{sp.split()[0][:12]}" for i,sp in zip(inv.isolate,inv.species)]
fig,ax=plt.subplots(figsize=(12, 0.7*len(cat_cols)+2))
Hn=(H/H.max().replace(0,1))  # column-normalised for colour; annotate raw counts
sns.heatmap(Hn.T, annot=H.T.astype(int), fmt="d", cmap="YlGnBu", cbar_kws={"label":"per-category max-normalised"},
            xticklabels=labels, yticklabels=[c.replace("SM:","") for c in cat_cols],
            linewidths=.5, linecolor="white", ax=ax, annot_kws={"size":8})
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title("Biocontrol-relevant gene inventory across UK yeast genomes\n"
             "(CAZyme lytic enzymes, MEROPS proteases, antiSMASH clusters)",
             fontsize=12, loc="left")
plt.setp(ax.get_xticklabels(), rotation=0, fontsize=8)
plt.savefig(f"{FIG}/fig11_biocontrol_genes.png", dpi=300, bbox_inches="tight"); plt.close()
print(f"\nWrote {DATA}/mining_gene_inventory.csv + {FIG}/fig11_biocontrol_genes.png ({len(inv)} genomes)")
