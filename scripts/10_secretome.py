#!/usr/bin/env python
"""
10_secretome.py
Restrict the biocontrol enzyme genotype to the SECRETED proteome (SignalP 5.0b),
so the gene counts are comparable to Gideon's plate assays, which only detect
enzymes the yeast actually exports.

For each mined genome (analysis/mining/results/<iso>/):
  <iso>.annotations.txt   -> per-transcript CAZyme family + MEROPS protease + Product
  <iso>.secreted_ids.txt  -> TranscriptIDs with a Sec/SPI signal peptide
The two are joined on TranscriptID (annotations col 2 == SignalP ID).

Outputs:
  analysis/data/secretome_inventory.csv      total vs secreted counts per class
  analysis/data/genotype_phenotype_secreted.csv
  analysis/figures/fig12_secretome.png       total vs secreted enzyme arsenal
"""
import os, re, glob, collections
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import spc, it

def pkey(sp):
    """Normalise mined species name to the shared palette key (figs 1-11)."""
    if sp.startswith("Metschnikowia"): return "Metschnikowia sp."
    if sp.startswith("Pichia"):        return "Pichia sp."
    return sp

MINE="analysis/mining"; DATA="analysis/data"; FIG="analysis/figures"
sheet=pd.read_csv(f"{MINE}/samplesheet.tsv", sep="\t")

# CAZy family -> biocontrol enzyme class (same mapping as 08_mining_summary.py)
CAZY_CLASS = {
 "Chitinase (GH18/19)":       ["GH18","GH19"],
 "b-1,3-glucanase (GH16/17/55/64/81/128)": ["GH16","GH17","GH55","GH64","GH81","GH128","GH152"],
 "b-1,6-glucanase (GH30)":    ["GH30"],
 "Cellulase (GH5/6/7/12/45)": ["GH5","GH6","GH7","GH12","GH45"],
 "Amylase (GH13/15)":         ["GH13","GH15"],
 "Mannanase (GH76)":          ["GH76"],
 "Pectinase (GH28/PL1/PL3/PL4)": ["GH28","PL1","PL3","PL4"],
}
FAM2CLASS={f:cls for cls,fl in CAZY_CLASS.items() for f in fl}

def classify(cazyme_cell):
    """Return set of enzyme classes a transcript's CAZyme annotation belongs to."""
    fams=re.findall(r"(GH\d+|GT\d+|CBM\d+|CE\d+|PL\d+|AA\d+)", str(cazyme_cell))
    return {FAM2CLASS[f] for f in fams if f in FAM2CLASS}

order=["RY14","RY13","TY26","TY30","UY1","TY6","TY4"]
rows=[]
for _,s in sheet.iterrows():
    iso=s.isolate; d=f"{MINE}/results/{iso}"
    ann=f"{d}/{iso}.annotations.txt"; secf=f"{d}/{iso}.secreted_ids.txt"
    if not (os.path.exists(ann) and os.path.exists(secf)):
        print("skip (missing):",iso); continue
    a=pd.read_csv(ann,sep="\t",dtype=str).fillna("")
    sec=set(l.strip() for l in open(secf) if l.strip())
    a["is_sec"]=a["TranscriptID"].isin(sec)
    rec=dict(isolate=iso, species=s.species.replace("_"," "),
             n_genes=len(a), n_secreted=int(a["is_sec"].sum()))
    # per enzyme class: total genes vs secreted genes
    for cls in CAZY_CLASS:
        m=a["CAZyme"].apply(lambda c: cls in classify(c))
        rec[f"{cls} | total"]=int(m.sum())
        rec[f"{cls} | secr"]=int((m & a["is_sec"]).sum())
    # MEROPS proteases
    prot=a["Protease"]!=""
    rec["MEROPS_protease | total"]=int(prot.sum())
    rec["MEROPS_protease | secr"]=int((prot & a["is_sec"]).sum())
    rows.append(rec)

inv=pd.DataFrame(rows)
inv["isolate"]=pd.Categorical(inv.isolate,[i for i in order if i in inv.isolate.values],ordered=True)
inv=inv.sort_values("isolate").reset_index(drop=True)
inv.to_csv(f"{DATA}/secretome_inventory.csv", index=False)
print(f"Wrote {DATA}/secretome_inventory.csv")
print(inv[["isolate","species","n_genes","n_secreted"]].to_string(index=False))
print("\nSecretome fraction: {:.1%}-{:.1%}".format(
    (inv.n_secreted/inv.n_genes).min(),(inv.n_secreted/inv.n_genes).max()))

# ---- FIGURE 12: secreted vs total enzyme arsenal ----
classes=list(CAZY_CLASS)+["MEROPS_protease"]
# drop classes that are all-zero total across genomes
classes=[c for c in classes if inv[f"{c} | total"].sum()>0]
labels=[f"{i}\n{sp.split()[0][:11]}" for i,sp in zip(inv.isolate,inv.species)]
ncol=len(inv); nrow=len(classes)
fig,axes=plt.subplots(nrow,1,figsize=(11,1.5*nrow+1),sharex=True)
if nrow==1: axes=[axes]
x=np.arange(ncol)
for ax,cls in zip(axes,classes):
    tot=inv[f"{cls} | total"].values; sec=inv[f"{cls} | secr"].values
    cols=[spc(pkey(sp)) for sp in inv.species]
    ax.bar(x,tot,color=cols,alpha=.30,width=.7,label="total")
    ax.bar(x,sec,color=cols,alpha=1.0,width=.7,label="secreted")
    for xi,(t,se) in enumerate(zip(tot,sec)):
        ax.text(xi,t+max(tot)*0.03+0.1,f"{se}/{t}",ha="center",va="bottom",fontsize=7)
    ax.set_ylabel(cls.split(" (")[0].replace("MEROPS_protease","Protease\n(MEROPS)"),
                  fontsize=8,rotation=0,ha="right",va="center")
    ax.set_ylim(0,max(tot.max(),1)*1.25)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
axes[0].set_title("Secreted (solid) vs total (faint) biocontrol enzyme arsenal\n"
                  "SignalP 5.0b Sec/SPI signal peptide — the fraction the plate assays can detect",
                  fontsize=11,loc="left")
axes[-1].set_xticks(x); axes[-1].set_xticklabels(labels,fontsize=8)
plt.tight_layout()
plt.savefig(f"{FIG}/fig12_secretome.png",dpi=300,bbox_inches="tight"); plt.close()
print(f"Wrote {FIG}/fig12_secretome.png")

# ---- SECRETED genotype vs phenotype (plate assays measure SECRETED enzymes) ----
RAW=f"{DATA}/raw_csv"
q=pd.read_csv(f"{RAW}/Lytic_enzyme_assay_results__Qualitative.csv",header=None,dtype=str).fillna("")
enz_names=["Cellulase","Amylase","Lipase","Chitinase","Protease","Esterase","Beta-glucosidase"]
qual={}
for _,r in q.iloc[2:].iterrows():
    iso=str(r[0]).strip()
    if not re.match(r"^(UY|TY|RY)\d+$",iso): continue
    calls={}
    for ei,en in enumerate(enz_names):
        c0=1+ei*3
        yes=sum(1 for k in range(3) if str(r[c0+k]).strip().lower()=="yes")
        calls[en]="+" if yes>=2 else "-"
    qual[iso]=calls
# map each assay to the secreted gene-count column
ASSAY2CLS={"Chitinase":"Chitinase (GH18/19)","Cellulase":"Cellulase (GH5/6/7/12/45)",
           "Amylase":"Amylase (GH13/15)","Protease":"MEROPS_protease",
           "Beta-glucosidase":"b-1,3-glucanase (GH16/17/55/64/81/128)"}
grows=[]
for _,g in inv.iterrows():
    iso=g.isolate; qc=qual.get(iso,{})
    rec=dict(isolate=iso,species=g.species)
    for assay,cls in ASSAY2CLS.items():
        rec[f"{assay}_secr"]=int(g[f"{cls} | secr"])
        rec[f"{assay}_total"]=int(g[f"{cls} | total"])
        rec[f"{assay}_assay"]=qc.get(assay,"n/a")
    grows.append(rec)
gp=pd.DataFrame(grows)
gp.to_csv(f"{DATA}/genotype_phenotype_secreted.csv",index=False)
print(f"\nWrote {DATA}/genotype_phenotype_secreted.csv")
print("\n=== SECRETED gene count (secr/total) | assay call ===")
for _,r in gp.iterrows():
    cells=" ".join(f"{a[:4].lower()} {r[a+'_secr']:>2}/{r[a+'_total']:<2}|{r[a+'_assay']}"
                   for a in ASSAY2CLS)
    print(f"  {r.isolate:5s} {r.species:26s} {cells}")
