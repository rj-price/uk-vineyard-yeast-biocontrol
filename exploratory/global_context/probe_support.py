import pandas as pd, numpy as np
from Bio import Phylo
GC = "analysis/global_context"
t = Phylo.read(f"{GC}/huvarum_phylo/supermatrix/huvarum_tree.contree", "newick")
names = pd.read_csv(f"{GC}/huvarum_names.csv"); kent = set(names[names.kind == "isolate"].name)

sup = [float(c.confidence) for c in t.get_nonterminals()
       if c.confidence is not None]
sup = np.array(sup)
print(f"internal nodes with support: {len(sup)}")
print(f"UFBoot: median {np.median(sup):.0f}, "
      f">=95: {(sup>=95).sum()}, 80-95: {((sup>=80)&(sup<95)).sum()}, <80: {(sup<80).sum()}")

# supermatrix length / informative sites
sm = f"{GC}/huvarum_phylo/supermatrix/SUPERMATRIX.fasta"
seqs = {}
name = None
for ln in open(sm):
    if ln.startswith(">"):
        name = ln[1:].strip(); seqs[name] = []
    else:
        seqs[name].append(ln.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
L = len(next(iter(seqs.values())))
arr = np.array([list(s) for s in seqs.values()])
var = sum(1 for j in range(L) if len(set(arr[:, j]) - {"-", "X"}) > 1)
print(f"supermatrix: {len(seqs)} taxa x {L} aa; variable columns ~{var} ({100*var/L:.1f}%)")
