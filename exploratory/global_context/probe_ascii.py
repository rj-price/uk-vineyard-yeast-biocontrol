import pandas as pd, io
from Bio import Phylo
GC = "analysis/global_context"
t = Phylo.read(f"{GC}/huvarum_phylo/supermatrix/huvarum_tree.contree", "newick")
names = pd.read_csv(f"{GC}/huvarum_names.csv"); kent = set(names[names.kind == "isolate"].name)
annot = pd.read_csv(f"{GC}/meta/huvarum_annot.csv"); cc = dict(zip(annot.acc, annot.country))
t.root_at_midpoint(); t.ladderize()
for cl in t.get_terminals():
    cl.name = ("*KENT_" + cl.name) if cl.name in kent else (cc.get(cl.name, "?")[:3] + "_" + cl.name)
buf = io.StringIO(); Phylo.draw_ascii(t, file=buf, column_width=64)
print(buf.getvalue())
