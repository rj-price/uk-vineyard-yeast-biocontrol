#!/usr/bin/env python
"""
11_comparative.py
Within-genus comparison of antagonistic vs non-antagonistic isolates (same species),
to find isolate-specific genomic differences that could explain the phenotype split.

Contrasts (dual-plate %MI vs B. cinerea):
  Rhodotorula babjevae : UY1 (Bc+ 34.5%) vs TY27/TY28/TY29 (Bc- 0%)   [iron/VOC ~equal]
  Hanseniaspora uvarum : RY14 (62%) > RY7 (36%) > RY25 (14%)          [tracks VOC+iron]

For every mined genome it parses (results/<iso>/):
  <iso>.annotations.txt   per-transcript CAZyme family, MEROPS protease, Product
  <iso>.secreted_ids.txt  SignalP Sec/SPI transcripts (secreted subset)
  antismash/*.json        secondary-metabolite cluster products
Outputs:
  data/comparative_<genus>.csv      family/cluster counts per isolate + assay
  figures/fig13_comparative.png     antagonist vs non differences, per genus
"""
import os, re, json, glob, collections
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(__file__))

MINE="analysis/mining"; DATA="analysis/data"; FIG="analysis/figures"
sheet=pd.read_csv(f"{MINE}/samplesheet_all.tsv", sep="\t")

# --- phenotype: dual/iron/VOC %MI vs B. cinerea ---
long=pd.read_csv(f"{DATA}/phenotypes_long.csv")
def mi(assay,iso,path="B. cinerea"):
    s=long[(long.assay==assay)&(long.isolate==iso)&(long.pathogen==path)]
    return round(s.pct_MI.mean(),1) if len(s) else np.nan

CAZY_CLASS = {
 "Chitinase (GH18/19)":["GH18","GH19"],
 "b-1,3-glucanase":["GH16","GH17","GH55","GH64","GH81","GH128","GH152"],
 "b-1,6-glucanase":["GH30"], "Cellulase":["GH5","GH6","GH7","GH12","GH45"],
 "Amylase":["GH13","GH15"], "Mannanase":["GH76"],
 "Pectinase":["GH28","PL1","PL3","PL4"],
}
FAM_RE=re.compile(r"(GH\d+|GT\d+|CBM\d+|CE\d+|PL\d+|AA\d+)")

def parse_genome(iso):
    d=f"{MINE}/results/{iso}"
    ann=f"{d}/{iso}.annotations.txt"; secf=f"{d}/{iso}.secreted_ids.txt"
    if not os.path.exists(ann): return None
    a=pd.read_csv(ann,sep="\t",dtype=str).fillna("")
    sec=set(l.strip() for l in open(secf)) if os.path.exists(secf) else set()
    a["is_sec"]=a["TranscriptID"].isin(sec)
    fam_tot=collections.Counter(); fam_sec=collections.Counter()
    for _,r in a.iterrows():
        fams=set(FAM_RE.findall(r["CAZyme"]))
        for f in fams:
            fam_tot[f]+=1
            if r["is_sec"]: fam_sec[f]+=1
    prot=a["Protease"]!=""
    # antiSMASH products
    prods=collections.Counter()
    js=glob.glob(f"{d}/antismash/*.json")
    if js:
        for rec in json.load(open(js[0])).get("records",[]):
            for ar in rec.get("areas",[]):
                for p in ar.get("products",[]): prods[p]+=1
    return dict(n_genes=len(a), n_secreted=int(a.is_sec.sum()),
                prot_tot=int(prot.sum()), prot_sec=int((prot&a.is_sec).sum()),
                fam_tot=fam_tot, fam_sec=fam_sec, prods=prods)

GROUPS={
 "Rhodotorula":dict(genus="Rhodotorula babjevae",
     order=["UY1","TY27","TY28","TY29"], antag=["UY1"]),
 "Hanseniaspora":dict(genus="Hanseniaspora uvarum",
     order=["RY14","RY7","RY25"], antag=["RY14"]),
}

def class_counts(g,which):
    src=g["fam_sec"] if which=="sec" else g["fam_tot"]
    return {cls:sum(src.get(f,0) for f in fl) for cls,fl in CAZY_CLASS.items()}

fig,axes=plt.subplots(1,2,figsize=(15,7))
for ax,(gname,cfg) in zip(axes,GROUPS.items()):
    isos=[i for i in cfg["order"] if parse_genome(i)]
    G={i:parse_genome(i) for i in isos}
    # build a tidy table: secreted lytic-enzyme classes + proteases + SM clusters
    rows=[]
    for i in isos:
        g=G[i]; cs=class_counts(g,"sec"); ct=class_counts(g,"tot")
        smb=collections.Counter()
        for p,c in g["prods"].items(): smb[p]+=c
        rec=dict(isolate=i, Bc_dual=mi("dual",i), Bc_iron=mi("iron",i), Bc_VOC=mi("VOC",i),
                 n_genes=g["n_genes"], n_secreted=g["n_secreted"],
                 protease_sec=g["prot_sec"], protease_tot=g["prot_tot"],
                 terpene=smb.get("terpene",0), NRPS=smb.get("NRPS",0)+smb.get("NRPS-like",0),
                 T1PKS=smb.get("T1PKS",0), siderophore=smb.get("NI-siderophore",0)+smb.get("NRPS-independent-siderophore",0))
        for cls in CAZY_CLASS:
            rec[f"{cls}|sec"]=cs[cls]; rec[f"{cls}|tot"]=ct[cls]
        rows.append(rec)
    df=pd.DataFrame(rows)
    df.to_csv(f"{DATA}/comparative_{gname}.csv", index=False)
    print(f"\n=== {gname} ({cfg['genus']}) — antagonist(s): {cfg['antag']} ===")
    cols=["isolate","Bc_dual","Bc_iron","Bc_VOC","n_secreted","protease_sec",
          "terpene","siderophore"]+[f"{c}|sec" for c in CAZY_CLASS]
    print(df[cols].to_string(index=False,na_rep="."))

    # figure: secreted lytic enzyme classes per isolate (the contact-mechanism candidates)
    plot_cls=[c for c in CAZY_CLASS if df[f"{c}|sec"].sum()>0]
    x=np.arange(len(plot_cls)); w=0.8/len(isos)
    for k,i in enumerate(isos):
        vals=[df.loc[df.isolate==i,f"{c}|sec"].values[0] for c in plot_cls]
        antag=i in cfg["antag"]
        bc=df.loc[df.isolate==i,"Bc_dual"].values[0]
        ax.bar(x+k*w, vals, w, label=f"{i} (Bc {bc:.0f}%)"+(" ★" if antag else ""),
               edgecolor="black" if antag else "none",
               linewidth=1.6 if antag else 0)
    ax.set_xticks(x+0.4-w/2); ax.set_xticklabels([c.split(" (")[0] for c in plot_cls],
                                                 rotation=30, ha="right", fontsize=9)
    ax.set_title(f"{cfg['genus']}\nsecreted lytic-enzyme arsenal (★ = Bc-antagonist)",
                 fontsize=11, loc="left")
    ax.set_ylabel("secreted genes"); ax.legend(fontsize=8, frameon=False)
    for s in ("top","right"): ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIG}/fig13_comparative.png", dpi=300, bbox_inches="tight"); plt.close()
print(f"\nWrote data/comparative_*.csv + {FIG}/fig13_comparative.png")
