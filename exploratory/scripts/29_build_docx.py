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
1: '**Fig. 1. Pathogen-specific antagonism across the vineyard yeast collection.** (A) Dual-culture mycelial inhibition (%MI) of the 25 isolates advanced from the qualitative pre-screen, against *Botrytis cinerea* Bc21010, *Fusarium oxysporum* f. sp. *lactucae* FoAJ516 and *Fusarium graminearum* Fg321, ordered by *B. cinerea* activity. Values are means of three biological replicates (n = 3); replicate values are given in Table S4. (B) UpSet plot of pathogen specificity, counting an isolate as active when its mean %MI was at least 5%. Most antagonists suppressed *B. cinerea* only; UY1, TY26 and TY30 were active against all three pathogens, and TY27, TY28 and TY29 against *F. oxysporum* only.',
2: '**Fig. 2. Phylogenomics and genome-based species identification.** (A) Maximum-likelihood phylogeny of the 55 sequenced isolates and reference and type-strain genomes, inferred from a concatenated supermatrix of 362 single-copy BUSCO loci (IQ-TREE, LG+F+G4, 1000 ultrafast bootstrap replicates). The tree is rooted on the Ascomycota/Basidiomycota split, which midpoint rooting recovers identically. Tips are coloured by genome-based species, and a superscript T marks reference assemblies from type material. (B) Whole-genome average nucleotide identity (ANI) of each isolate against the type and reference genomes of its own genus, for the five genera assigned by ANI. The 95% species threshold is marked on the colour scale, and cells marked <76 fell below the fastANI reporting limit. *Vishniacozyma* and *Filobasidium* were assigned from D1/D2 rDNA sequence (Table S1); for *Vishniacozyma* the only available reference genomes are from non-type strains, and TY7 shares 85.4% ANI with the *V. victoriae* reference. The cross-genus tree resolves relationships within species less well than the *H. uvarum* tree in Fig. 5B. The cryptic second *Sporobolomyces* species is shown in full in Fig. S3.',
3: '**Fig. 3. Complementary antagonistic mechanisms.** (A) Volatile-mediated inhibition (%MI) in the sealed split-plate assay. (B) Iron dose response, showing %MI at 5 mg l⁻¹ (filled) and 20 mg l⁻¹ (open) FeCl₃. The assay had no unsupplemented arm, so inference rests on the 5 to 20 mg l⁻¹ comparison only. All eight *Fusarium* pairings declined (sign test, *P* = 0.008, n = 8), whereas the eleven *B. cinerea* pairings did not (sign test, *P* = 0.23, n = 11). (C) Relative extracellular lytic-enzyme activity (clear-zone/colony ratio) across the enzyme panel; a ratio of 1.0 means no halo. Because the ratio scales halo to colony size, it penalises fast-growing isolates. Values are means of three biological replicates (n = 3) for the isolates advanced to the mechanistic assays.',
4: '**Fig. 4. Genomic basis of biocontrol traits.** (A) Polysaccharide-degrading CAZyme genes by class in the 13 mined antagonist genomes (dbCAN, genes supported by at least two methods). (B) Total and SignalP-predicted secreted proteases (MEROPS). (C) Biosynthetic gene clusters (antiSMASH) by type, including the pulcherrimin cyclodipeptide synthase (CDPS) found only in *Metschnikowia*. (D) Within-genus contrasts in *Rhodotorula* aff. *babjevae* (UY1, TY27–TY29) and *Hanseniaspora uvarum* (RY14, RY7 and RY25). The left panel shows dual-culture %MI against each pathogen. The right panel shows total CAZymes, total and secreted proteases and biosynthetic gene clusters for the same isolates, shaded by deviation from the genus mean; a difference of one gene is within annotation noise and is left unshaded. *B. cinerea* activity varied within both genera, whereas gene counts changed little.',
5: '**Fig. 5. Global genomic context of the two lead biocontrol lineages.** (A) All-versus-all ANI across the *Rhodotorula babjevae*/*graminis* complex. The Kent isolates form a distinct third species-level cluster with public genomes from China, Australia and an Italian culture collection that were deposited as *R. graminis* or *R. babjevae*. The cluster falls below the 95% threshold against both the *R. babjevae* (~90.3%) and *R. graminis* (~88.8%) type strains. (B) BUSCO-supermatrix maximum-likelihood phylogeny of the 22 Kent *H. uvarum* isolates and 16 global conspecifics, midpoint-rooted; re-deposits of the same strain are shown once. Global tips are annotated by country and isolation source. Kent tips are filled if the isolate suppressed *B. cinerea* in the dual-culture screen and open if it was not active in the screen. The Kent isolates fall into several clades of mostly Kent isolates; all nine members of the lowest Kent clade (100% UFBoot) are *B. cinerea* antagonists. Scale bar, substitutions per site.',
6: '**Fig. 6. Candidate prioritisation.** Highest-ranked isolates, coloured by species, with their pathogen %MI and mechanism profile (volatile %MI, mean %MI on FeCl₃-supplemented medium and enzyme breadth). Cells are shaded as a fraction of the column maximum. The composite score is the sum of *z*-scores for cross-pathogen %MI (weighted ×2), *B. cinerea* %MI, volatile %MI, mean %MI on FeCl₃-supplemented medium and enzyme breadth (Methods). RY14 ranks first under 69% of 242 alternative weightings and in 99% of replicate bootstrap resamples; below the top three, ranks are not stable. The score orders candidates for follow-up within this dataset and does not predict efficacy.',
}
FIGFILE = {1:"Fig1_screen.png",2:"Fig2_phylogenomics.png",3:"Fig3_mechanisms.png",
           4:"Fig4_genomic.png",5:"Fig5_global.png",6:"Fig6_prioritisation.png"}

def figblock(n):
    return f"\n{img(FIGFILE[n])}\n\n{LEG[n]}\n"

# ---- manuscript with figures inserted before the next subsection ----
man = strip_comments(open(f"{MAN}/manuscript_draft.md").read())
ANCHOR = {
1: "### Whole-genome phylogeny and species identity",
2: "### Complementary antagonistic mechanisms",
3: "### Genomic basis of biocontrol traits",
4: "### Kent lineages are cosmopolitan and partly mislabelled in public databases",
5: "### Candidate prioritisation",
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
# the assembly-rescue figure was dropped in revision round 3
for n, f in {1:"FigS1_clonal_ani.png", 2:"FigS2_biocontrol_global.png",
             3:"FigS3_sporobolomyces_ani.png"}.items():
    supp = supp.replace(f"**Fig. S{n}.", f"{img(f)}\n\n**Fig. S{n}.", 1)

# ---- assemble ----
cover = strip_comments(open(f"{MAN}/cover_letter.md").read())
doc = cover + PB + man + PB + "# Supplementary Material\n\n" + \
      supp.split("# Supplementary Material", 1)[1]
open(f"{MAN}/submission.md", "w").write(doc)
print("wrote submission.md  chars=", len(doc))

# ---- pandoc -> docx, then add continuous line numbering (Access Microbiology requires it) ----
import subprocess, zipfile, shutil, tempfile, re as _re

DOCX = f"{MAN}/submission.docx"
PANDOC = "/mnt/apps/users/jnprice/conda/envs/earlgrey/bin/pandoc"
subprocess.run([PANDOC, f"{MAN}/submission.md",
                "-o", DOCX, "--resource-path", BASE], check=True)


def add_line_numbers(path):
    """Insert <w:lnNumType> into every section so the .docx carries continuous line
    numbering. Word only renders line numbers when the section properties declare them,
    and pandoc has no option for this, so it is patched into the built file."""
    ln = '<w:lnNumType w:countBy="1" w:restart="continuous" w:distance="360"/>'
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        z.extractall(tmp)
    dp = os.path.join(tmp, "word", "document.xml")
    xml = open(dp, encoding="utf-8").read()
    if "lnNumType" not in xml:
        # Within a sectPr the schema orders lnNumType after footnotePr/endnotePr/type/
        # pgSz/pgMar, so insert it just before </w:sectPr>, which satisfies that order
        # whichever of the optional earlier elements pandoc emitted.
        xml, n = _re.subn(r"</w:sectPr>", ln + "</w:sectPr>", xml)
        print(f"line numbering added to {n} section(s)")
    open(dp, "w", encoding="utf-8").write(xml)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for nm in names:
            z.write(os.path.join(tmp, nm), nm)
    shutil.rmtree(tmp)


add_line_numbers(DOCX)
print("wrote submission.docx")
