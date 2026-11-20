import pandas as pd, os
from Bio import Phylo
GC = "analysis/global_context"
t = Phylo.read(f"{GC}/huvarum_phylo/supermatrix/huvarum_tree.contree", "newick")

names = pd.read_csv(f"{GC}/huvarum_names.csv")
annot = pd.read_csv(f"{GC}/meta/huvarum_annot.csv")
kent = set(names[names.kind == "isolate"].name)
acc_country = dict(zip(annot.acc, annot.country))
acc_src = dict(zip(annot.acc, annot.source))


def lab(n):
    if n in kent:
        return f"KENT:{n}"
    return f"{acc_country.get(n,'?')}|{acc_src.get(n,'?')[:12]}:{n}"


# midpoint root for readability
t.root_at_midpoint()
t.ladderize()
print("=== ladderised tip order (top→bottom) ===")
for lf in t.get_terminals():
    print("  ", lf.name in kent and ">>" or "  ", lf.name, "|", lab(lf.name))

# are all Kent isolates monophyletic?
kent_tips = [lf for lf in t.get_terminals() if lf.name in kent]
mrca = t.common_ancestor(kent_tips)
under = [lf.name for lf in mrca.get_terminals()]
non_kent_under = [x for x in under if x not in kent]
print(f"\nKent MRCA subtends {len(under)} tips; non-Kent inside = {non_kent_under}")
print("Kent monophyletic:", len(non_kent_under) == 0)
