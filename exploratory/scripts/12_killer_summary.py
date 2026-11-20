#!/usr/bin/env python
"""
12_killer_summary.py
Summarise the killer-toxin tblastn (analysis/mining/killer_toxin/killer_tblastn.tsv):
curated fungal killer-toxin panel vs all 60 assemblies. Classify each hit as a
genuine toxin match vs a conserved-domain cross-hit, and cross-reference antagonism.
-> analysis/data/killer_toxin_hits.csv
"""
import pandas as pd, numpy as np
DATA="analysis/data"; KT="analysis/mining/killer_toxin/killer_tblastn.tsv"
cols=['qseqid','qlen','sseqid','pident','length','mismatch','gapopen','evalue','bitscore','qcovs']
d=pd.read_csv(KT,sep='\t',names=cols)
d['isolate']=d.sseqid.str.replace('__.*','',regex=True)
d['toxin']=d.qseqid.str.split('|').str[2]
m=pd.read_csv(f"{DATA}/isolate_master.csv")[['isolate','species']]
long=pd.read_csv(f"{DATA}/phenotypes_long.csv")
def bc(iso):
    s=long[(long.assay=='dual')&(long.isolate==iso)&(long.pathogen=='B. cinerea')]
    return round(s.pct_MI.mean(),1) if len(s) else np.nan

# genuine toxin = decent identity AND covers a real fraction of the (mature) toxin
d['call']=np.where((d.pident>=30)&(d.qcovs>=40),'genuine',
          np.where(d.qcovs<30,'domain cross-hit','weak'))
d=d.merge(m,on='isolate',how='left')
d['Bc_dual']=d.isolate.map(bc)
d=d.sort_values(['call','bitscore'],ascending=[True,False])
keep=['isolate','species','Bc_dual','toxin','pident','length','qcovs','evalue','bitscore','call']
d[keep].to_csv(f"{DATA}/killer_toxin_hits.csv",index=False)

print("=== GENUINE killer-toxin hits (pident>=30, qcovs>=40) ===")
g=d[d.call=='genuine']
print(g[keep].to_string(index=False) if len(g) else "(none)")
print("\n=== zymocin-alpha (KTXA) = endochitinase(GH18) domain cross-hit — NOT a toxin ===")
z=d[d.toxin=='KTXA_KLULA']
print(f"  {z.isolate.nunique()} isolates, pident {z.pident.min():.0f}-{z.pident.max():.0f}%, "
      f"qcovs {z.qcovs.min()}-{z.qcovs.max()}% (partial = chitinase domain only)")
print("\nToxins with NO genuine hit anywhere:")
allt={'M11P_SCVM1','TOXK_MILFA','KHS1_YEASX','TOXK_CYBMR','KP4T_UMV4','KP6T_UMV6',
      'TOXK_CYBSA','KTXG_KLULA'}
print("  ",sorted(allt-set(g.toxin)))
print(f"\nWrote {DATA}/killer_toxin_hits.csv")
