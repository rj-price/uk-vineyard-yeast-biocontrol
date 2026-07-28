#!/usr/bin/env python
"""
01_build_master.py
Build the master isolate key linking:
  genome code (CUY/CTY/RY) <-> assay code (UY/TY/RY) <-> sequencing id (Y1_/Y2_/Y3_)
  <-> species <-> assembly stats / BUSCO <-> phenotype assays.

Outputs (analysis/data/):
  isolate_master.csv     one row per isolate: mapping, species, assembly, BUSCO
  phenotype_matrix.csv   isolate x trait wide table (numeric) for PCA/clustering/heatmap
  phenotypes_long.csv    tidy long form of all %MI assays with replicates
"""
import re, os
import numpy as np
import pandas as pd

RAW = "analysis/data/raw_csv"
OUT = "analysis/data"
os.makedirs(OUT, exist_ok=True)

def rd(name, **kw):
    return pd.read_csv(os.path.join(RAW, name), dtype=str, keep_default_na=False, **kw)

# ----------------------------------------------------------------------------
# 1. Genomics sheet -> mapping + assembly + BUSCO  (rows 3-64 in original)
# ----------------------------------------------------------------------------
g = rd("Yeast_Assemblies_and_Analysis__Sheet1.csv", header=None)
cols = ["genome_code","seq_id","morph_genome","its1_blast",
        "size_inclN","size_noN","seq_num","mean_size","median_size",
        "longest","shortest","gc","n50","l50","n90","gap",
        "busco_complete","busco_single","busco_dup","busco_frag","busco_missing"]
g = g.iloc[2:64].copy()          # data rows only (skip 2 headers, stop before NOTE block)
g.columns = cols + [f"x{i}" for i in range(len(g.columns)-len(cols))]
g = g[g["genome_code"].str.match(r"^(CUY|CTY|RY)\d+$", na=False)].copy()

# assay code = genome code without leading culture "C"
def to_assay(code):
    if code.startswith("CUY"): return "UY"+code[3:]
    if code.startswith("CTY"): return "TY"+code[3:]
    return code                              # RY unchanged
g["isolate"] = g["genome_code"].apply(to_assay)

num = ["size_inclN","size_noN","seq_num","mean_size","median_size","longest",
       "shortest","gc","n50","l50","n90","gap","busco_complete","busco_single",
       "busco_dup","busco_frag","busco_missing"]
for c in num:
    g[c] = pd.to_numeric(g[c], errors="coerce")

# ----------------------------------------------------------------------------
# 2. Species clean-up (genome ITS1 BLAST + BUSCO-tree clustering notes)
# ----------------------------------------------------------------------------
def clean_species(s):
    s = s.strip()
    if s.startswith("Rhodosporidium babjevae"): return "Rhodotorula babjevae"
    if s.startswith("Sporobolomyces roseus"):   return "Sporobolomyces roseus"
    if s.startswith("Vishniacozyma carnescens"):return "Vishniacozyma carnescens"
    if s.startswith("Naganishia"):              return "Naganishia vaughanmartiniae"
    if s.startswith("Filobasidium"):            return "Filobasidium stepposum"
    if s.startswith("Hanseniaspora uvarum"):    return "Hanseniaspora uvarum"
    if "Hannaella" in s:                        return "Hannaella sp."
    if "Metschnikowia" in s:                    return "Metschnikowia sp."
    if "Pichia" in s:                           return "Pichia sp."
    if "clusters with H. uvarum" in s:          return "Hanseniaspora uvarum"
    if "clusters with Y3_4" in s:               return "Pichia sp."
    if "clusters with Y2_6" in s:               return "Sporobolomyces sp."
    if s == "NO ITS1":                          return "Metschnikowia sp."   # Y2_9 red-halo, sister to Y2_13
    return s or "unresolved"
g["species"] = g["its1_blast"].apply(clean_species)

# ----------------------------------------------------------------------------
# 3. Assay-side identity / morphology (Combined Data has best per-isolate ITS)
# ----------------------------------------------------------------------------
cd = rd("Combined_Data__Sheet1.csv", header=None)
cd = cd.iloc[4:].copy()
cd.columns = (["isolate","colour","texture","cell_shape","its_assay",
               "Fg321","FoAJ516","Bc21010","Bc21024","BcCHF","VdRH1020",
               "VdRH1022","VdRH1024","Fg_ave","Fg_sd","Fo_ave","Fo_sd",
               "Bc_ave","Bc_sd","barcode","organism"]
              + [f"y{i}" for i in range(len(cd.columns)-21)])
cd = cd[cd["isolate"].str.match(r"^(UY|TY|RY)\d+$", na=False)].copy()
assay_id = cd.set_index("isolate")[["colour","texture","cell_shape","its_assay","barcode"]]

# General classification morphology (covers all isolates, incl. non-sequenced)
gc = rd("Classification_of_yeast__isolates__General_classification.csv", header=None)
gc = gc.iloc[2:].copy()
gc.columns = ["isolate","col_colour","col_texture","cell_shape"] + \
             [f"z{i}" for i in range(len(gc.columns)-4)]
gc = gc[gc["isolate"].str.match(r"^(UY|TY|RY)\d+$", na=False)].set_index("isolate")

# ----------------------------------------------------------------------------
# 4. Master table
# ----------------------------------------------------------------------------
m = g[["isolate","genome_code","seq_id","species","its1_blast","morph_genome",
       "size_inclN","gc","n50","l50","seq_num","longest",
       "busco_complete","busco_dup","busco_frag","busco_missing"]].copy()
m = m.merge(assay_id, left_on="isolate", right_index=True, how="left")
m = m.merge(gc[["col_colour","col_texture","cell_shape"]].rename(
              columns={"col_colour":"colour_gc","col_texture":"texture_gc",
                       "cell_shape":"cell_shape_gc"}),
            left_on="isolate", right_index=True, how="left")
# prefer assay ITS species where it is more precise than genome
m["species_assayITS"] = m["its_assay"].replace("", np.nan)
m = m.sort_values("seq_id").reset_index(drop=True)
m.to_csv(os.path.join(OUT,"isolate_master.csv"), index=False)

# ----------------------------------------------------------------------------
# 5. Phenotype assays -> long + wide
# ----------------------------------------------------------------------------
def read_mi(fname, assay, pathogen, condition=""):
    df = rd(fname, header=None)
    df = df[df[0].str.match(r"^(UY|TY|RY)\d+$", na=False)].iloc[:, :6].copy()
    df.columns = ["isolate","R1","R2","R3","mean_growth","pct_MI"]
    for c in ["R1","R2","R3","mean_growth","pct_MI"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["assay"]=assay; df["pathogen"]=pathogen; df["condition"]=condition
    return df

long_parts = [
    read_mi("Dual_plate_assays__B_cinerea.csv","dual","B. cinerea"),
    read_mi("Dual_plate_assays__F_graminearum.csv","dual","F. graminearum"),
    read_mi("Dual_plate_assays__F_oxysporum.csv","dual","F. oxysporum"),
    read_mi("VOC_Assay__B_cinerea.csv","VOC","B. cinerea"),
    read_mi("VOC_Assay__F_graminearum.csv","VOC","F. graminearum"),
    read_mi("VOC_Assay__F_oxysporum.csv","VOC","F. oxysporum"),
    read_mi("Iron_competition_assay__B_cinerea_+_5_mgL.csv","iron","B. cinerea","5 mg/L Fe"),
    read_mi("Iron_competition_assay__B_cinerea_+_20_mgL.csv","iron","B. cinerea","20 mg/L Fe"),
    read_mi("Iron_competition_assay__F_graminearum_+_5_mgL.csv","iron","F. graminearum","5 mg/L Fe"),
    read_mi("Iron_competition_assay__F_graminearum_+_20_mgL.csv","iron","F. graminearum","20 mg/L Fe"),
    read_mi("Iron_competition_assay__F_oxysporum_+_5_mgL.csv","iron","F. oxysporum","5 mg/L Fe"),
    read_mi("Iron_competition_assay__F_oxysoprum_+_20_mgL.csv","iron","F. oxysporum","20 mg/L Fe"),
]
long = pd.concat(long_parts, ignore_index=True)
long.to_csv(os.path.join(OUT,"phenotypes_long.csv"), index=False)

# ----------------------------------------------------------------------------
# 6. Enzyme quantitative (relative activity) -> per enzyme
# ----------------------------------------------------------------------------
eq = rd("Lytic_enzyme_assay_results__Quantitative.csv", header=None)
enzymes = ["Cellulase","Amylase","Lipase","Protease","Esterase","Chitinase","Beta-glucosidase"]
cur=None; rows=[]
for _,r in eq.iterrows():
    tag = str(r[1]).strip()
    for e in enzymes:
        if tag.startswith(e): cur=e
    if re.match(r"^(UY|TY|RY)\d+$", str(r[0])):
        rows.append((str(r[0]).strip(), cur, pd.to_numeric(r[9], errors="coerce")))
enz = pd.DataFrame(rows, columns=["isolate","enzyme","rel_activity"])
enz_wide = enz.pivot_table(index="isolate", columns="enzyme", values="rel_activity")
enz_wide.columns = [f"enz_{c}" for c in enz_wide.columns]

# ----------------------------------------------------------------------------
# 7. Wide phenotype matrix for multivariate analysis
# ----------------------------------------------------------------------------
def mi_pivot(assay, tag):
    sub = long[long["assay"]==assay].copy()
    sub["key"] = tag + "_" + sub["pathogen"].str.replace(". ","",regex=False).str.replace(" ","")
    if assay=="iron":
        sub["key"] = sub["key"] + "_" + sub["condition"].str.replace(" ","").str.replace("/","")
    return sub.pivot_table(index="isolate", columns="key", values="pct_MI")

wide = pd.concat([mi_pivot("dual","MI"), mi_pivot("VOC","VOC"),
                  mi_pivot("iron","Fe"), enz_wide], axis=1)
wide = wide.merge(m.set_index("isolate")[["species","seq_id","genome_code"]],
                  left_index=True, right_index=True, how="left")
wide.index.name="isolate"
wide.to_csv(os.path.join(OUT,"phenotype_matrix.csv"))

# ----------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------
print("=== isolate_master.csv ===", m.shape)
print(m["species"].value_counts().to_string())
print("\nSequenced isolates with phenotype data:",
      wide.dropna(subset=["MI_Bcinerea"]).shape[0] if "MI_Bcinerea" in wide else "n/a")
print("\n=== phenotype_matrix.csv columns ===")
print([c for c in wide.columns if c not in ("species","seq_id","genome_code")])
print("\nRows in phenotype_matrix:", wide.shape[0])
print("\nSpecies of the phenotyped (dual-plate) isolates:")
ph = wide[wide.index.isin(long[long.assay=='dual'].isolate.unique())]
print(ph[["species"]].dropna().sort_index().to_string())
