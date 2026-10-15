#!/usr/bin/env python
"""
36_accessory.py
Accessory gene presence/absence within H. uvarum (22 Kent isolates) and Rhodotorula
(8 Kent isolates), from the miniprot projection in allelic/08_accessory.sh
(reviewer point 1).

Call per pan-gene x genome, best miniprot hit:
  absent    = no hit covering >= 50% of the protein at >= 90% identity
  disrupted = hit present but with a frameshift or in-frame stop
  intact    = otherwise
A gene is "variable" if its state is not the same in every genome of the set.

Tests
  H. uvarum: each variable gene vs B. cinerea activity in the dual-culture screen
  (Fisher exact, BH-corrected), and whether it simply marks the active 9-isolate clade.
  Rhodotorula: genes whose state differs inside a near-clonal group
  (A = UY1, TY27, UY4; B = TY28, UY2, UY3; fastANI ~99.99%).

-> analysis/data/accessory_calls_<genus>.csv   (pan-gene x genome states)
-> analysis/data/accessory_summary.csv
"""
import glob, os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis"
ACC = f"{PROJ}/allelic/accessory"
CLADE = {"RY3", "RY4", "RY5", "RY8", "RY19", "RY20", "RY22", "RY25", "RY26"}
GROUPS = {"A": ["UY1", "TY27", "UY4"], "B": ["TY28", "UY2", "UY3"]}


def calls(genus):
    pan = [l[1:].split()[0] for l in open(f"{ACC}/{genus}/pan.faa") if l.startswith(">")]
    cols = {}
    for f in sorted(glob.glob(f"{ACC}/{genus}/paf/*.paf")):
        g = os.path.basename(f)[:-4]
        best = {}
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            q, ql, qs, qe = x[0], int(x[1]), int(x[2]), int(x[3])
            ident = int(x[9]) / max(int(x[10]), 1)
            tags = dict(t.split(":", 2)[::2] for t in x[12:] if t.count(":") >= 2)
            cov = (qe - qs) / ql
            sc = int(tags.get("AS", 0))
            if q not in best or sc > best[q][0]:
                bad = int(tags.get("fs", 0)) > 0 or int(tags.get("st", 0)) > 0
                best[q] = (sc, cov, ident, bad)
        st = {}
        for q in pan:
            b = best.get(q)
            if b is None or b[1] < 0.5 or b[2] < 0.9:
                st[q] = "absent"
            elif b[3]:
                st[q] = "disrupted"
            else:
                st[q] = "intact"
        cols[g] = st
    d = pd.DataFrame(cols).loc[pan]
    d.to_csv(f"{PROJ}/data/accessory_calls_{genus}.csv")
    return d


rows = []
pre = pd.read_csv(f"{PROJ}/data/prelim_screen.csv").set_index("isolate")

# ---- H. uvarum
h = calls("Hanseniaspora")
intact = (h == "intact").astype(int)
var = intact[(intact.sum(axis=1) > 0) & (intact.sum(axis=1) < intact.shape[1])]
act = (pre.loc[h.columns, "GI_Bc_ave"].fillna(0) > 0).astype(int)
clade = pd.Series([int(c in CLADE) for c in h.columns], index=h.columns)
ps, perfect_act, perfect_clade = [], 0, 0
for gid, r in var.iterrows():
    t = pd.crosstab(r, act).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    ps.append(fisher_exact(t.values).pvalue)
    perfect_act += int((r == act).all() or (r == 1 - act).all())
    perfect_clade += int((r == clade).all() or (r == 1 - clade).all())
q = multipletests(ps, method="fdr_bh")[1] if ps else []
print(f"H. uvarum: {len(h)} pan-genes, {len(var)} variable (intact in some, not all, of {h.shape[1]})")
print(f"  min BH q vs B. cinerea activity = {min(q):.3f}; genes q<0.05: {(q < 0.05).sum()}")
print(f"  genes perfectly tracking activity: {perfect_act}; perfectly marking the 9-isolate clade: {perfect_clade}")
rows.append(dict(set="H. uvarum (22)", pan_genes=len(h), variable=len(var),
                 q_min=round(min(q), 4), n_q05=int((q < 0.05).sum()),
                 perfect_activity=perfect_act, perfect_clade=perfect_clade))
top = pd.DataFrame({"gene": var.index, "p": ps, "q": q}).sort_values("p").head(15)
top = top.join(var.loc[top.gene].reset_index(drop=True).set_axis(top.index), how="left")
print(top[["gene", "p", "q"]].to_string(index=False))

# ---- Rhodotorula
r = calls("Rhodotorula")
for g, mem in GROUPS.items():
    sub = r[mem]
    diff = sub[sub.nunique(axis=1) > 1]
    hard = diff[(diff == "intact").any(axis=1) & (diff == "absent").any(axis=1)]
    print(f"Rhodotorula group {g} {mem}: {len(diff)} genes differ in state, "
          f"{len(hard)} intact in one member and absent in another")
    if len(hard):
        print(hard.to_string())
    rows.append(dict(set=f"Rhodotorula group {g} ({', '.join(mem)})", pan_genes=len(r),
                     variable=len(diff), intact_vs_absent=len(hard)))
pd.DataFrame(rows).to_csv(f"{PROJ}/data/accessory_summary.csv", index=False)
