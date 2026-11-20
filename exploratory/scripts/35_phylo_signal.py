#!/usr/bin/env python
"""
35_phylo_signal.py
Phylogenetic signal of B. cinerea antagonism within Kent H. uvarum (reviewer point 1).

Trait: B. cinerea effect in Gideon's qualitative competition pre-screen (any of
Bc21010, Bc21024, BcCHF), which covers all 22 sequenced Kent H. uvarum; the
quantitative dual-culture screen covers only the 13 advanced isolates.
Test: Fitch parsimony steps on the BUSCO supermatrix tree pruned to the Kent tips,
against 10,000 tip-label permutations. Fewer steps than the null = signal.

-> analysis/data/phylo_signal.csv
"""
import numpy as np, pandas as pd
from Bio import Phylo

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis"
TREE = f"{PROJ}/global_context/huvarum_phylo/supermatrix/huvarum_tree.contree"
rng = np.random.default_rng(1)
NPERM = 10000

pre = pd.read_csv(f"{PROJ}/data/prelim_screen.csv").set_index("isolate")
huv = [l.strip() for l in open(f"{PROJ}/global_context/huvarum_isolates.txt") if l.strip()]

tree = Phylo.read(TREE, "newick")
for tip in [t for t in tree.get_terminals() if t.name not in huv]:
    tree.prune(tip)
tips = [t.name for t in tree.get_terminals()]
assert sorted(tips) == sorted(huv), (tips, huv)


def fitch(tree, state):
    """Fitch parsimony steps for a binary trait; state: tip name -> 0/1."""
    steps = 0

    def rec(clade):
        nonlocal steps
        if clade.is_terminal():
            return {state[clade.name]}
        sets = [rec(c) for c in clade.clades]
        cur = sets[0]
        for s in sets[1:]:
            inter = cur & s
            if inter:
                cur = inter
            else:
                cur = cur | s
                steps += 1
        return cur
    rec(tree.root)
    return steps


rows = []
traits = {
    "B. cinerea (pre-screen, any strain)": (pre.loc[tips, ["Bc21010", "Bc21024", "BcCHF"]].max(axis=1) > 0).astype(int),
    "B. cinerea (dual culture, %MI > 0; untested = 0)": (pre.loc[tips, "GI_Bc_ave"].fillna(0) > 0).astype(int),
}
for name, s in traits.items():
    st = s.to_dict()
    obs = fitch(tree, st)
    vals = np.array(list(st.values()))
    null = np.empty(NPERM, int)
    for i in range(NPERM):
        null[i] = fitch(tree, dict(zip(tips, rng.permutation(vals))))
    p = (np.sum(null <= obs) + 1) / (NPERM + 1)
    rows.append(dict(trait=name, n=len(tips), n_active=int(vals.sum()), steps_obs=obs,
                     steps_null_mean=null.mean().round(2), steps_null_min=null.min(), p_perm=p))
    print(name, "active", vals.sum(), "/", len(tips), "steps", obs, "null mean", null.mean().round(2), "p", round(p, 5))

# the one well-supported Kent-only clade (100% UFBoot in Fig 5B); Fisher test inside vs outside
from scipy.stats import fisher_exact
CLADE = {"RY3", "RY4", "RY5", "RY8", "RY19", "RY20", "RY22", "RY25", "RY26"}
for name, s in traits.items():
    a = int(s[list(CLADE)].sum()); b = int(s.drop(list(CLADE)).sum()); nb = len(tips) - len(CLADE)
    p = fisher_exact([[a, len(CLADE) - a], [b, nb - b]]).pvalue
    rows.append(dict(trait=name + " | clade vs rest (Fisher)", n=len(tips), n_active=a + b,
                     clade_active=f"{a}/{len(CLADE)}", rest_active=f"{b}/{nb}", p_fisher=round(p, 5)))
    print(name, "clade", a, "/", len(CLADE), "rest", b, "/", nb, "Fisher p", round(p, 5))

print("tip order:", ", ".join(f"{t}:{traits[list(traits)[0]][t]}" for t in tips))
pd.DataFrame(rows).to_csv(f"{PROJ}/data/phylo_signal.csv", index=False)
