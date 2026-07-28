#!/usr/bin/env python
"""
27_supp_tables.py -> analysis/manuscript/supplementary/Table_S1..S10.csv
Assembles submission-ready supplementary tables from the project data files.
Kent isolate set = 55 sequenced (60 genome codes minus the 5 cryptic Sporobolomyces
TY8/TY9/TY12/TY17/TY24 dropped from the paper).
"""
import os, re, glob
import numpy as np, pandas as pd

D = "analysis/data"
OUT = "analysis/manuscript/supplementary"
os.makedirs(OUT, exist_ok=True)
DROP = {"TY8", "TY9", "TY12", "TY17", "TY24"}

SRC = {"U": "Chardonnay (untreated)", "T": "Chardonnay (fungicide-treated)", "R": "Pinot Noir"}
def berry(iso):
    return SRC.get(iso[0], "")

def genus(sp):
    return str(sp).split()[0] if pd.notna(sp) else ""

master = pd.read_csv(f"{D}/isolate_master.csv")
final = pd.read_csv(f"{D}/species_id_final.csv")
seq = master[~master.isolate.isin(DROP)].copy()          # 55 sequenced
fin = final[~final.isolate.isin(DROP)].copy()

def order(df):
    df = df.copy()
    df["_p"] = df.isolate.str[0].map({"U": 0, "T": 1, "R": 2})
    df["_n"] = df.isolate.str.extract(r"(\d+)").astype(int)
    return df.sort_values(["_p", "_n"]).drop(columns=["_p", "_n"])

# ---------- Table S1 : isolates + genome-based identification ----------
s1 = order(fin.merge(seq[["isolate", "genome_code", "seq_id"]], on="isolate", how="left"))
s1["berry_source"] = s1.isolate.map(berry)
s1["id_method"] = np.where(s1.basis.astype(str).str.startswith("ANI"),
                           "Whole-genome ANI", "rDNA D1/D2 (LSU)")
s1 = s1[["isolate", "genome_code", "seq_id", "berry_source", "prior_label",
         "species_genome", "ani_best_ref", "ani", "id_method", "confidence",
         "genus_changed", "species_changed"]]
s1.columns = ["Isolate", "Genome code", "Sequencing ID", "Berry source", "Prior (ITS) identification",
              "Genome-based species", "Closest reference", "ANI (%)", "ID method", "Confidence",
              "Genus revised", "Species revised"]
s1["ANI (%)"] = s1["ANI (%)"].round(2)
s1.to_csv(f"{OUT}/Table_S1_isolates.csv", index=False)

# ---------- Table S2 : assembly & BUSCO statistics ----------
s2 = order(seq.copy())
s2["Assembler"] = np.where(s2.species_genome.map(genus).isin(["Hanseniaspora", "Pichia"]),
                           "Platanus-allee", "SPAdes")
s2["Assembly (Mb)"] = (s2.size_inclN / 1e6).round(2)
s2["N50 (kb)"] = (s2.n50 / 1e3).round(1)
s2 = s2[["isolate", "species_genome", "Assembler", "Assembly (Mb)", "seq_num", "N50 (kb)",
         "l50", "gc", "busco_complete", "busco_dup", "busco_frag", "busco_missing"]]
s2.columns = ["Isolate", "Species", "Assembler", "Assembly (Mb)", "Contigs", "N50 (kb)", "L50",
              "GC (%)", "BUSCO complete (%)", "BUSCO duplicated (%)", "BUSCO fragmented (%)",
              "BUSCO missing (%)"]
for c in ["GC (%)", "BUSCO complete (%)", "BUSCO duplicated (%)", "BUSCO fragmented (%)", "BUSCO missing (%)"]:
    s2[c] = s2[c].round(1)
s2.to_csv(f"{OUT}/Table_S2_assembly_busco.csv", index=False)

# ---------- Table S3 : reference / type-strain genomes used for ANI & phylogeny ----------
rows = []
for f in sorted(glob.glob("analysis/species_id/refs/*.fna")):
    b = os.path.basename(f)[:-4]
    m = re.match(r"(.+?)__(GC[AF]_\d+\.\d+)", b)
    if m:
        rows.append({"Species / strain": m.group(1).replace("_", " "), "Assembly accession": m.group(2)})
s3 = pd.DataFrame(rows).sort_values("Species / strain")
s3.to_csv(f"{OUT}/Table_S3_reference_genomes.csv", index=False)

# ---------- Table S4 : antagonism screen (quantitative %MI for the 25 antagonists) ----------
ps = pd.read_csv(f"{D}/prelim_screen.csv")
ps = ps[ps[["GI_Bc_ave", "GI_Fo_ave", "GI_Fg_ave"]].notna().any(axis=1)].copy()
ps = ps.merge(fin[["isolate", "species_genome"]], on="isolate", how="left")
ps["B. cinerea %MI (mean of 3 isolates ± SD)"] = ps.apply(
    lambda r: f"{r.GI_Bc_ave:.1f} ± {r.GI_Bc_sd:.1f}" if pd.notna(r.GI_Bc_ave) else "", axis=1)
ps["F. oxysporum %MI (± SD)"] = ps.apply(
    lambda r: f"{r.GI_Fo_ave:.1f} ± {r.GI_Fo_sd:.1f}" if pd.notna(r.GI_Fo_ave) else "", axis=1)
ps["F. graminearum %MI (± SD)"] = ps.apply(
    lambda r: f"{r.GI_Fg_ave:.1f} ± {r.GI_Fg_sd:.1f}" if pd.notna(r.GI_Fg_ave) else "", axis=1)
s4 = order(ps)[["isolate", "species_genome",
                "B. cinerea %MI (mean of 3 isolates ± SD)", "F. oxysporum %MI (± SD)",
                "F. graminearum %MI (± SD)"]]
s4.columns = ["Isolate", "Genome-based species", "B. cinerea %MI (mean of Bc21010/Bc21024/BcCHF ± SD)",
              "F. oxysporum FoAJ516 %MI (± SD)", "F. graminearum Fg321 %MI (± SD)"]
s4.to_csv(f"{OUT}/Table_S4_screen.csv", index=False)

# ---------- Table S5 : mechanistic assays ----------
pm = pd.read_csv(f"{D}/phenotype_matrix.csv")
enzc = [c for c in pm.columns if c.startswith("enz_")]
s5 = order(pm)[["isolate", "species",
                "VOC_Bcinerea", "VOC_Fgraminearum", "VOC_Foxysporum",
                "Fe_Bcinerea_5mgLFe", "Fe_Bcinerea_20mgLFe",
                "Fe_Fgraminearum_5mgLFe", "Fe_Fgraminearum_20mgLFe",
                "Fe_Foxysporum_5mgLFe", "Fe_Foxysporum_20mgLFe"] + enzc]
s5.columns = ["Isolate", "Species",
              "VOC %MI B. cinerea", "VOC %MI F. graminearum", "VOC %MI F. oxysporum",
              "%MI B. cinerea +5 mg/L Fe", "%MI B. cinerea +20 mg/L Fe",
              "%MI F. graminearum +5 mg/L Fe", "%MI F. graminearum +20 mg/L Fe",
              "%MI F. oxysporum +5 mg/L Fe", "%MI F. oxysporum +20 mg/L Fe"] + \
             [c.replace("enz_", "Enzyme (REA): ") for c in enzc]
s5 = s5.round(2)
s5.to_csv(f"{OUT}/Table_S5_mechanistic.csv", index=False)

# ---------- Table S6 : biocontrol gene inventory of mined genomes ----------
mi = pd.read_csv(f"{D}/mining_gene_inventory.csv")
sec = pd.read_csv(f"{D}/secretome_inventory.csv")
s6 = mi.merge(sec[["isolate", "n_secreted", "MEROPS_protease | secr"]], on="isolate", how="left")
s6 = s6.rename(columns={"n_genes": "Predicted genes", "n_secreted": "Secreted proteins (SignalP)",
                        "MEROPS_protease": "MEROPS proteases (total)",
                        "MEROPS_protease | secr": "MEROPS proteases (secreted)",
                        "n_SM_clusters": "antiSMASH SM clusters", "n_cazyme": "CAZymes (total)"})
s6.to_csv(f"{OUT}/Table_S6_gene_inventory.csv", index=False)

# ---------- Table S7 : killer-toxin tblastn hits ----------
kt = pd.read_csv(f"{D}/killer_toxin_hits.csv")
kt = kt[~kt.isolate.isin(DROP)]
s7 = order(kt)[["isolate", "species", "toxin", "pident", "length", "qcovs", "evalue", "call", "Bc_dual"]]
s7.columns = ["Isolate", "Species", "Killer toxin (query)", "Identity (%)", "Aln length",
              "Query cov (%)", "E-value", "Call", "B. cinerea antagonist"]
s7.to_csv(f"{OUT}/Table_S7_killer_toxin.csv", index=False)

# ---------- Table S8 : clonal-pair (UY1/TY27) protein-altering variants ----------
cv = pd.read_csv(f"{D}/clonal_variants_ipr.csv")
s8 = cv[["contig", "pos", "ref", "alt", "kind", "effect", "aa_change", "transcript",
         "product", "cazyme", "protease", "secreted", "pfam", "ipr_desc", "go", "ANY_regulatory"]].copy()
s8.columns = ["Contig", "Position", "Ref", "Alt", "Type", "Effect", "AA change", "Transcript",
              "Product", "CAZyme", "Protease", "Secreted", "PFAM", "InterPro description",
              "GO terms", "Regulatory"]
s8.to_csv(f"{OUT}/Table_S8_clonal_variants.csv", index=False)

# ---------- Table S9 : global comparator genomes ----------
def load_annot(f, sp_col):
    a = pd.read_csv(f)
    a = a.rename(columns={"acc": "Assembly accession", "sp": "Species (as deposited)",
                          "strain": "Strain", "level": "Assembly level", "country": "Country",
                          "source": "Isolation source category", "iso_src": "Isolation source (raw)",
                          "date": "Year"})
    keep = ["Assembly accession", "Species (as deposited)", "Strain", "Country",
            "Isolation source (raw)", "Year", "Assembly level"]
    return a[[c for c in keep if c in a.columns]]

h = load_annot(f"analysis/global_context/meta/huvarum_annot.csv", "sp")
r = load_annot(f"analysis/global_context/meta/rhodotorula_annot.csv", "sp")
h.insert(0, "Set", "H. uvarum (global)")
r.insert(0, "Set", "Rhodotorula genus")
s9 = pd.concat([h, r], ignore_index=True)
s9.to_csv(f"{OUT}/Table_S9_global_comparators.csv", index=False)

# ---------- Table S10 : biocontrol gene content, Kent vs global ----------
bg = pd.read_csv(f"{D}/biocontrol_global_inventory.csv")
bg.to_csv(f"{OUT}/Table_S10_biocontrol_kent_vs_global.csv", index=False)

# ---------- summary ----------
for f in sorted(glob.glob(f"{OUT}/Table_S*.csv")):
    df = pd.read_csv(f)
    print(f"{os.path.basename(f):42s} rows={len(df):4d}  cols={len(df.columns)}")
