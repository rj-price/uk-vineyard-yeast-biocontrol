#!/usr/bin/env python
"""
29_build_docx.py -> analysis/manuscript/submission.md  (then pandoc -> submission.docx)
Assembles one submission document: cover letter -> manuscript with Fig 1-6 + legends
integrated -> supplementary (Fig S1-S3 + legends + Tables S1-S10 + legends).
"""
import os, re, glob
import pandas as pd
from PIL import Image

BASE = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
MAN = f"{BASE}/analysis/manuscript"
FIGDIR = f"{BASE}/analysis/figures/manuscript"
SUPP = f"{MAN}/supplementary"

def strip_comments(t):
    return re.sub(r"<!--.*?-->", "", t, flags=re.S).strip() + "\n"

def img(fname, maxw=6.3, maxh=8.0):
    p = f"{FIGDIR}/{fname}"
    w, h = Image.open(p).size
    width = min(maxw, round(maxh * w / h, 2))
    return f"![]({p}){{width={width}in}}"

PB = '\n```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n'

# ---- main-figure legends ----
LEG = {
1: "**Fig. 1. Pathogen-specific antagonism across the vineyard yeast collection.** (A) Quantitative dual-culture mycelial inhibition (%MI) for the 25 isolates advanced from the antagonism screen, against *Botrytis cinerea* (Bc21010), *Fusarium oxysporum* (FoAJ516) and *Fusarium graminearum* (Fg321); rows ordered by *B. cinerea* activity, with a genome-based species colour bar. (B) UpSet plot of pathogen specificity: most antagonists suppress *B. cinerea* only, and only three are active across pathogen genera. *Verticillium dahliae* was included in preliminary screening but omitted (minimal activity).",
2: "**Fig. 2. Phylogenomics and genome-based species identification.** (A) Maximum-likelihood phylogeny inferred from a concatenated single-copy BUSCO supermatrix (IQ-TREE, LG+F+G4, 1000 ultrafast bootstrap replicates) of the 55 sequenced isolates with reference/type-strain genomes and a *Fusarium venenatum* outgroup; tips coloured by genome-based species. (B) Average nucleotide identity (ANI) heatmap (isolates × references) with the 95 % species threshold marked, resolving the sub-threshold *Rhodotorula* lineage, the *Pichia kluyveri* assignment, the *Metschnikowia pulcherrima* complex, and the cryptic second *Sporobolomyces* species.",
3: "**Fig. 3. Complementary antagonistic mechanisms.** (A) Volatile-mediated inhibition (%MI). (B) Iron competition: %MI without and with FeCl₃ supplementation. (C) Relative extracellular lytic-enzyme activity (clear-zone/colony ratio) across the enzyme panel. Values are for the isolates advanced to mechanistic assays.",
4: "**Fig. 4. Genomic basis of biocontrol traits.** (A) Cell-wall-degrading CAZyme inventory (counts per class) across the 13 mined antagonist genomes, grouped by species. (B) Protease arsenal: total versus SignalP-secreted counts. (C) antiSMASH secondary-metabolite clusters by species, including the pulcherrimin cyclodipeptide synthase (CDPS) in *Metschnikowia*. (D) KT395-type killer-toxin presence across the *Pichia kluyveri* isolates. (E) Within-genus phenotype for the two conspecific contrasts (*Rhodotorula* cf. *babjevae* UY1/TY27–TY29; *Hanseniaspora uvarum* RY14/RY7/RY25): per-pathogen %MI, with *B. cinerea* activity varying while *Fusarium* activity is conserved. (F) Matched gene content (secretome, secondary-metabolite clusters, secreted enzymes) shown as deviation from the genus mean, flat across the phenotypically divergent isolates.",
5: "**Fig. 5. Global genomic context of the two lead biocontrol lineages.** (A) All-versus-all ANI across the *Rhodotorula babjevae*/*graminis* complex: the Kent isolates form a distinct third species-level cluster (with mislabelled public genomes from China and Australia), separate from both the *R. babjevae* and *R. graminis* type strains. (B) BUSCO-supermatrix phylogeny of the Kent *Hanseniaspora uvarum* isolates with global conspecifics, tips annotated by country and source: the Kent isolates are interspersed with worldwide grape/wine strains, indicating a cosmopolitan species with no distinct UK lineage.",
6: "**Fig. 6. From lab to fruit.** (A) Exploratory detached grape-berry assay, one representative berry per condition at 10 days post-inoculation (water only; yeast only, *H. uvarum* RY14; *B. cinerea* only; *B. cinerea* + yeast). Water- and yeast-only berries remained infection-free; both *B. cinerea*-inoculated conditions developed grey mould, with no clear protection. The assay was qualitative (no disease scoring). (B) Genome-characterised candidate-prioritisation dashboard: top isolates ranked by a composite mechanism score alongside their pathogen %MI and mechanism profile (volatiles, iron competition, secreted-enzyme breadth), coloured by species.",
}
FIGFILE = {1:"Fig1_screen.png",2:"Fig2_phylogenomics.png",3:"Fig3_mechanisms.png",
           4:"Fig4_genomic.png",5:"Fig5_global.png",6:"Fig6_labtofruit.png"}

def figblock(n):
    return f"\n{img(FIGFILE[n])}\n\n{LEG[n]}\n"

# ---- manuscript with figures inserted before the next subsection ----
man = strip_comments(open(f"{MAN}/manuscript_draft.md").read())
ANCHOR = {
1: "### Whole-genome phylogeny and species identity (Fig 2)",
2: "### Complementary antagonistic mechanisms (Fig 3)",
3: "### Genomic basis of biocontrol traits (Fig 4)",
4: "### The clonal pair narrows, but does not resolve, the determinant",
5: "### From lab to fruit: the translational gap and candidate prioritisation (Fig 6)",
6: "## Discussion",
}
for n, anc in ANCHOR.items():
    assert anc in man, f"anchor not found: {anc}"
    man = man.replace(anc, figblock(n) + "\n" + anc, 1)

# ---- supplementary: figures inserted before their legends; tables are NOT embedded ----
# Per author decision, the supplementary tables are supplied as separate CSV files
# (Table_S1..Table_S10.csv) rather than embedded in the document; only their legends
# appear here. This keeps the submission document from ballooning with long tables.
supp = strip_comments(open(f"{SUPP}/supplementary_materials.md").read())
for n, f in {1:"FigS1_clonal_ani.png", 2:"FigS2_biocontrol_global.png",
             3:"FigS3_assembly_rescue.png", 4:"FigS4_sporobolomyces_ani.png"}.items():
    supp = supp.replace(f"**Fig. S{n}.", f"{img(f)}\n\n**Fig. S{n}.", 1)

# ---- assemble ----
cover = strip_comments(open(f"{MAN}/cover_letter.md").read())
doc = cover + PB + man + PB + "# Supplementary Material\n\n" + \
      supp.split("# Supplementary Material", 1)[1]
open(f"{MAN}/submission.md", "w").write(doc)
print("wrote submission.md  chars=", len(doc))
