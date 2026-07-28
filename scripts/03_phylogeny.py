#!/usr/bin/env python
"""
03_phylogeny.py
Render the BUSCO supermatrix ML tree (IQ-TREE) with tips relabelled to
isolate code + species, Fusarium venenatum (C42A/B) as outgroup root,
tips coloured by species. -> analysis/figures/fig6_phylogeny.png
Also writes a relabelled Newick for downstream use.
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import SP_ORDER as sp_all, SP_COL as sp_col

DATA="analysis/data"; FIG="analysis/figures"
TREE="genomics/busco_phylogenomics/out/supermatrix/SUPERMATRIX.phylip.treefile"

master = pd.read_csv(f"{DATA}/isolate_master.csv")
seq2info = {r.seq_id:(r.isolate, r.species) for _,r in master.iterrows()}

tree = Phylo.read(TREE,"newick")

def token(name):
    m = re.match(r"^(C42[AB]|Y\d+_\d+)", name or "")
    return m.group(1) if m else name

# relabel + colour
tip_species={}
for t in tree.get_terminals():
    tok = token(t.name)
    if tok in ("C42A","C42B"):
        iso, sp = tok, "Fusarium venenatum (outgroup)"
    else:
        iso, sp = seq2info.get(tok,(tok,"unresolved"))
    tip_species[iso if iso else tok] = sp
    t.name = f"{iso}  {sp if 'outgroup' not in sp else 'F. venenatum (outgroup)'}"
    t._sp = sp

# root on Fusarium outgroup
og = [t for t in tree.get_terminals() if t.name.startswith(("C42A","C42B"))]
try:
    tree.root_with_outgroup(*og)
except Exception:
    pass
tree.ladderize()

# collapse near-zero support display; draw
fig, ax = plt.subplots(figsize=(10,13))
def label_colour(clade):
    return "black"
def do_label(clade):
    return clade.name if clade.is_terminal() else ""
Phylo.draw(tree, axes=ax, do_show=False, label_func=do_label,
           show_confidence=False, label_colors=None)
# recolour tip label text by species
for txt in ax.texts:
    s = txt.get_text().strip()
    for sp,c in sp_col.items():
        key = "F. venenatum (outgroup)" if "outgroup" in sp else sp
        if s.endswith(key):
            txt.set_color(c); txt.set_fontsize(7.5); break
ax.set_title("BUSCO supermatrix phylogeny of UK yeast isolates\n"
             "(50 taxa, 1,759 aa, IQ-TREE VT+G4; rooted on F. venenatum)",
             fontsize=12, loc="left")
ax.set_xlabel("substitutions/site"); ax.set_ylabel("")
ax.set_xlim(-0.02, 1.55)                 # extra right margin for tip labels
handles=[Patch(facecolor=sp_col[s],
               label=(f"$\\it{{{s.replace(' ',chr(92)+' ')}}}$" if "sp." not in s and "outgroup" not in s else s))
         for s in sp_all]
handles.append(Patch(facecolor="#000000",label="F. venenatum (outgroup)"))
ax.legend(handles=handles, frameon=False, fontsize=8.5, title="Species",
          title_fontsize=9.5, loc="center left", bbox_to_anchor=(0.02, 0.42))
plt.savefig(f"{FIG}/fig6_phylogeny.png", dpi=300, bbox_inches="tight"); plt.close()

Phylo.write(tree, f"{DATA}/relabelled_tree.nwk", "newick")
print("Tips:", tree.count_terminals())
print("Species on tree:", pd.Series(list(tip_species.values())).value_counts().to_string())
print("Wrote fig6_phylogeny.png and relabelled_tree.nwk")
