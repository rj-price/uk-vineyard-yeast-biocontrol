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
1: "**Fig. 1. Pathogen-specific antagonism across the vineyard yeast collection.** (A) Qualitative effect/no-effect screen of the 56 isolates against five pathogen strains (three *Botrytis cinerea*: Bc21010, Bc21024, BcCHF; *Fusarium oxysporum* FoAJ516; *Fusarium graminearum* Fg321), rows grouped by genome-based species. (B) Quantitative dual-culture mycelial inhibition (%MI) for the 25 isolates advanced from the screen, against *B. cinerea* (mean of three strains), *F. oxysporum* and *F. graminearum*. (C) UpSet plot of pathogen specificity: most antagonists suppress *B. cinerea* only, and only three are active across pathogen genera. *Verticillium dahliae* was included in preliminary screening but omitted (minimal activity).",
2: "**Fig. 2. Phylogenomics and genome-based species identification.** (A) Maximum-likelihood phylogeny inferred from a concatenated single-copy BUSCO supermatrix (IQ-TREE, LG+F+G4, 1000 ultrafast bootstrap replicates) of the 55 sequenced isolates with reference/type-strain genomes and a *Fusarium venenatum* outgroup; tips coloured by genome-based species. (B) Average nucleotide identity (ANI) heatmap (isolates × references) with the 95 % species threshold marked, resolving the sub-threshold *Rhodotorula* lineage, the *Pichia kluyveri* assignment, the *Metschnikowia pulcherrima* complex, and the cryptic second *Sporobolomyces* species.",
3: "**Fig. 3. Complementary antagonistic mechanisms.** (A) Volatile-mediated inhibition (%MI). (B) Iron competition: %MI without and with FeCl₃ supplementation. (C) Relative extracellular lytic-enzyme activity (clear-zone/colony ratio) across the enzyme panel. Values are for the isolates advanced to mechanistic assays.",
4: "**Fig. 4. Genomic basis of biocontrol traits.** (A) Biocontrol-gene inventory (cell-wall-degrading CAZyme classes, MEROPS proteases, secondary-metabolite clusters) across the seven representative mined genomes. (B) Total versus SignalP-secreted enzyme and protease counts. (C) antiSMASH secondary-metabolite clusters by species, including the pulcherrimin cyclodipeptide synthase (CDPS) in *Metschnikowia*. (D) KT395-type killer-toxin presence across the *Pichia kluyveri* isolates.",
5: "**Fig. 5. *Botrytis*-specific antagonism without a gene-content basis.** (A) Secreted-enzyme genotype versus qualitative assay concordance. (B) The two within-species contrasts (*Rhodotorula* cf. *babjevae* UY1 vs TY27–TY29; *Hanseniaspora uvarum* RY14/RY7/RY25): per-pathogen %MI (replicate points) above matched biocontrol gene-content bars, showing gene content flat while *B. cinerea* activity varies. (C) Within-group ANI, highlighting the near-clonal UY1/TY27 pair (99.99 %). (D) The 188 protein-altering variants distinguishing the clonal pair, by functional class.",
6: "**Fig. 6. From lab to fruit.** (A) Exploratory detached grape-berry assay, one representative berry per condition at 10 days post-inoculation (water only; yeast only, *H. uvarum* RY14; *B. cinerea* only; *B. cinerea* + yeast). Water- and yeast-only berries remained infection-free; both *B. cinerea*-inoculated conditions developed grey mould, with no clear protection. The assay was qualitative (no disease scoring). (B) Genome-characterised candidate-prioritisation dashboard: top isolates ranked by a composite mechanism score alongside their pathogen %MI and mechanism profile (volatiles, iron competition, secreted-enzyme breadth), coloured by species.",
}
FIGFILE = {1:"Fig1_screen.png",2:"Fig2_phylogenomics.png",3:"Fig3_mechanisms.png",
           4:"Fig4_genomic.png",5:"Fig5_within.png",6:"Fig6_labtofruit.png"}

def figblock(n):
    return f"\n{img(FIGFILE[n])}\n\n{LEG[n]}\n"

# ---- manuscript with figures inserted before the next subsection ----
man = strip_comments(open(f"{MAN}/manuscript_draft.md").read())
ANCHOR = {
1: "### Whole-genome phylogeny and species identity (Fig 2)",
2: "### Complementary antagonistic mechanisms (Fig 3)",
3: "### Genomic basis of biocontrol traits (Fig 4)",
4: "### *Botrytis*-specific antagonism is not encoded by gene content (Fig 5)",
5: "### A conserved-pathway, regulatory explanation in *H. uvarum*",
6: "## Discussion",
}
for n, anc in ANCHOR.items():
    assert anc in man, f"anchor not found: {anc}"
    man = man.replace(anc, figblock(n) + "\n" + anc, 1)

# ---- supplementary: figures before legends, tables after legends ----
supp = strip_comments(open(f"{SUPP}/supplementary_materials.md").read())
for n, f in {1:"FigS1_global_context.png",2:"FigS2_biocontrol_global.png",3:"FigS3_assembly_rescue.png"}.items():
    supp = supp.replace(f"**Fig. S{n}.", f"{img(f)}\n\n**Fig. S{n}.", 1)

def md_table(csv):
    df = pd.read_csv(csv).fillna("").astype(str)
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]).replace("|", "\\|").replace("\n", " ") for c in cols) + " |")
    return "\n".join(out)

tblfile = {int(re.search(r"Table_S(\d+)_", f).group(1)): f
           for f in glob.glob(f"{SUPP}/Table_S*.csv")}
lines, out = supp.split("\n"), []
for ln in lines:
    out.append(ln)
    m = re.match(r"\*\*Table S(\d+)\.", ln)
    if m:
        n = int(m.group(1))
        out.append("")
        out.append(md_table(tblfile[n]))
        out.append("")
supp = "\n".join(out)

# ---- assemble ----
cover = strip_comments(open(f"{MAN}/cover_letter.md").read())
doc = cover + PB + man + PB + "# Supplementary Material\n\n" + \
      supp.split("# Supplementary Material", 1)[1]
open(f"{MAN}/submission.md", "w").write(doc)
print("wrote submission.md  chars=", len(doc))
