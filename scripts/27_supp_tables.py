#!/usr/bin/env python
"""
27_supp_tables.py -> analysis/manuscript/supplementary/Table_S1..S10.csv
Assembles submission-ready supplementary tables from the project data files.
Kent isolate set = the 56 screened isolates, 55 of them sequenced.
"""
import os, re, glob, json
import numpy as np, pandas as pd

D = "analysis/data"
OUT = "analysis/manuscript/supplementary"
os.makedirs(OUT, exist_ok=True)
SCREENED = set(pd.read_csv(f"{D}/prelim_screen.csv").isolate)

SRC = {"U": "Chardonnay (untreated)", "T": "Chardonnay (fungicide-treated)", "R": "Pinot Noir"}
def berry(iso):
    return SRC.get(iso[0], "")

def genus(sp):
    return str(sp).split()[0] if pd.notna(sp) else ""

master = pd.read_csv(f"{D}/isolate_master.csv")
final = pd.read_csv(f"{D}/species_id_final.csv")
seq = master[master.isolate.isin(SCREENED)].copy()          # 55 sequenced
fin = final[final.isolate.isin(SCREENED)].copy()

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
s1 = s1[["isolate", "genome_code", "seq_id", "berry_source",
         "species_genome", "ani_best_ref", "ani", "id_method", "confidence"]]
s1.columns = ["Isolate", "Genome code", "Sequencing ID", "Berry source",
              "Genome-based species", "Closest reference", "ANI (%)", "ID method", "Confidence"]
s1["ANI (%)"] = s1["ANI (%)"].round(2)
s1.to_csv(f"{OUT}/Table_S1_isolates.csv", index=False)
# Every later table takes its species label from S1, so names stay consistent.
SP = dict(zip(s1.Isolate, s1["Genome-based species"]))
def sp(iso):
    return SP.get(iso, "not sequenced")

# ---------- Table S2 : assembly & BUSCO statistics ----------
s2 = order(seq.copy())
s2["Assembler"] = np.where(s2.species_genome.map(genus).isin(["Hanseniaspora", "Pichia"]),
                           "Platanus-allee", "SPAdes")
# isolate_master holds the ORIGINAL SPAdes statistics for every isolate. The Hanseniaspora
# and Pichia genomes that were re-assembled with Platanus-allee are the ones actually
# annotated, mined and deposited, so their rows must carry the post-Platanus statistics
# (analysis/data/platanus_comparison.csv), not the SPAdes ones they were reported with.
_pl = pd.read_csv(f"{D}/platanus_comparison.csv").set_index("isolate")
_m = s2.isolate.isin(_pl.index) & (s2.Assembler == "Platanus-allee")
for col, src in [("size_inclN", "new_size"), ("n50", "new_n50"), ("seq_num", "new_contigs"),
                 ("busco_complete", "new_busco"), ("busco_dup", "new_D")]:
    s2.loc[_m, col] = s2.loc[_m, "isolate"].map(_pl[src]).values
# Full BUSCO breakdown for the re-assembled genomes comes from their own short summaries.
import glob, re as _re
_b = {}
for _f in glob.glob("analysis/reassembly/platanus/*/short_summary.*.txt"):
    _sid = os.path.basename(os.path.dirname(_f))
    _mm = _re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%",
                     open(_f).read())
    if _mm:
        _b[_sid] = dict(zip(("C", "S", "D", "F", "M"), map(float, _mm.groups())))
_bd = pd.DataFrame(_b).T
for col, src in [("busco_complete", "C"), ("busco_dup", "D"),
                 ("busco_frag", "F"), ("busco_missing", "M")]:
    s2.loc[_m, col] = s2.loc[_m, "seq_id"].map(_bd[src]).values
# L50 is not carried in the comparison file; recompute it from the deposited assembly.
def _l50(path):
    lens = []
    cur = 0
    for ln in open(path):
        if ln.startswith(">"):
            if cur:
                lens.append(cur)
            cur = 0
        else:
            cur += len(ln.strip())
    if cur:
        lens.append(cur)
    lens.sort(reverse=True)
    tot = sum(lens); run = 0
    for i, L in enumerate(lens):
        run += L
        if run >= tot / 2:
            return i + 1
    return np.nan
s2.loc[_m, "l50"] = [
    _l50(f"analysis/final_assemblies/ascomycetes_platanus/{sid}.fasta")
    for sid in s2.loc[_m, "seq_id"]]
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
# Species, strain and type status come from NCBI (datasets summary genome, saved in
# species_id/ref_metadata.tsv), not from the download file names, some of which are stale.
meta = pd.read_csv("analysis/species_id/ref_metadata.tsv", sep="\t", header=None,
                   names=["acc", "organism", "strain", "type", "level"], keep_default_na=False)
meta = meta.set_index("acc")
rows = []
for f in sorted(glob.glob("analysis/species_id/refs/*.fna")):
    m = re.match(r"(.+?)__(GC[AF]_\d+\.\d+)", os.path.basename(f)[:-4])
    if m:
        r_ = meta.loc[m.group(2)]
        if "heterotypic" in r_.type:
            tm = "yes (heterotypic synonym)"
        elif r_.type.startswith("assembly from type"):
            tm = "yes"
        else:
            tm = "no"
        rows.append({"Species (NCBI)": " ".join(r_.organism.split()[:2]), "Strain": r_.strain,
                     "Type material": tm, "Assembly accession": m.group(2),
                     "Assembly level": r_.level, "Used for": "species ID and phylogeny"})
# M. pulcherrima-clade genomes compared with TY26/TY30 (species_id/02_metsch_clade.sh)
SYN = {"NRRL Y-27328": "M. fructicola", "277": "M. fructicola", "CBS 10809T": "M. andauensis"}
mani = pd.read_csv("analysis/species_id/metsch_clade/ani.tsv", sep="\t", header=None,
                   names=["q", "r", "ani", "m", "t"])
mani["acc"] = mani.r.str.extract(r"/(GC[AF]_\d+\.\d+)/")[0]
mani["iso"] = mani.q.str.extract(r"/(TY\d+)\.fasta")[0]
seen = {r["Assembly accession"] for r in rows}
for line in open("analysis/species_id/metsch_clade/meta.jsonl"):
    j = json.loads(line)
    acc = j["accession"]
    a = mani[mani.acc == acc].set_index("iso").ani
    if a.empty:
        continue
    strain = (j["organism"].get("infraspecific_names") or {}).get("strain", "")
    lab = (j.get("type_material") or {}).get("type_label", "")
    tm = (f"yes (heterotypic synonym; type of {SYN[strain]})" if lab == "TYPE_MATERIAL_SYN"
          else "yes" if lab == "TYPE_MATERIAL" else "no")
    org = j["organism"]["organism_name"].split()
    sp_ = " ".join(org[:3] if org[1] == "aff." else org[:2])
    note = f"pulcherrima-clade ANI (TY26 {a.get('TY26', np.nan):.1f}%, TY30 {a.get('TY30', np.nan):.1f}%)"
    if acc in seen:
        for r in rows:
            if r["Assembly accession"] == acc:
                r["Used for"] += "; " + note
        continue
    rows.append({"Species (NCBI)": sp_, "Strain": strain, "Type material": tm,
                 "Assembly accession": acc, "Assembly level": j["assembly_info"]["assembly_level"],
                 "Used for": note})
s3 = pd.DataFrame(rows).sort_values(["Species (NCBI)", "Strain"])
s3.to_csv(f"{OUT}/Table_S3_reference_genomes.csv", index=False)

# ---------- Table S4 : quantitative dual-culture screen, individual replicates ----------
# Raw colony diameters (mm) of the three biological replicates per yeast-pathogen pair,
# from the dual-plate sheets. Controls reached the plate edge (85 mm) in every replicate,
# so per-replicate %MI = (85 - D) / 85 x 100. B. cinerea is strain Bc21010 throughout.
DC = 85.0
PATH = [("B. cinerea Bc21010", "B_cinerea", "MI_Bcinerea"),
        ("F. oxysporum f. sp. lactucae FoAJ516", "F_oxysporum", "MI_Foxysporum"),
        ("F. graminearum Fg321", "F_graminearum", "MI_Fgraminearum")]
_pm = pd.read_csv(f"{D}/phenotype_matrix.csv").set_index("isolate")
rows = []
for lab, f, pmcol in PATH:
    raw = pd.read_csv(f"{D}/raw_csv/Dual_plate_assays__{f}.csv").iloc[:, :4]
    raw.columns = ["isolate", "R1", "R2", "R3"]
    ctrl = raw[raw.isolate == "Control"][["R1", "R2", "R3"]].values
    assert (ctrl == DC).all(), f"{f}: control diameter is not {DC} mm"
    for _, r_ in raw[raw.isolate != "Control"].iterrows():
        # replicate means must reproduce the %MI used in the figures
        mi = (DC - r_[["R1", "R2", "R3"]].mean()) / DC * 100
        assert abs(mi - _pm.loc[r_.isolate, pmcol]) < 0.01, (r_.isolate, pmcol)
        for k in ("R1", "R2", "R3"):
            rows.append({"isolate": r_.isolate, "Pathogen": lab, "Replicate": int(k[1]),
                         "Colony diameter (mm)": r_[k],
                         "%MI": round((DC - r_[k]) / DC * 100, 1)})
s4 = pd.DataFrame(rows)
s4["_g"] = s4.isolate.str[0].map({"U": 0, "T": 1, "R": 2})
s4["_i"] = s4.isolate.str.extract(r"(\d+)", expand=False).astype(int)
s4["_k"] = s4.Pathogen.map({l: i for i, (l, _, _) in enumerate(PATH)})
s4 = s4.sort_values(["_g", "_i", "_k", "Replicate"]).drop(columns=["_g", "_i", "_k"])
s4.insert(1, "Genome-based species", s4.isolate.map(sp))
s4 = s4.rename(columns={"isolate": "Isolate"})
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
              "%MI B. cinerea, 5 mg/l FeCl3", "%MI B. cinerea, 20 mg/l FeCl3",
              "%MI F. graminearum, 5 mg/l FeCl3", "%MI F. graminearum, 20 mg/l FeCl3",
              "%MI F. oxysporum, 5 mg/l FeCl3", "%MI F. oxysporum, 20 mg/l FeCl3"] + \
             [c.replace("enz_", "Enzyme (REA): ") for c in enzc]
s5 = s5[s5.iloc[:, 2:].notna().any(axis=1)].copy()      # only isolates actually assayed
s5["Species"] = s5.Isolate.map(sp)
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
# fungal-wall CAZymes (chitinase, b-1,3/b-1,6-glucanase, mannanase) on the same dbCAN
# >=2-tool call as the class columns, total and with a SignalP signal peptide
WALL = {"GH18", "GH19", "GH16", "GH17", "GH55", "GH64", "GH81", "GH128", "GH152", "GH30", "GH76"}
def wall_counts(iso):
    d = f"analysis/mining/results/{iso}"
    ov = pd.read_csv(f"{d}/dbcan/overview.tsv", sep="\t")
    ov.columns = [c.strip() for c in ov.columns]
    ov = ov[pd.to_numeric(ov["#ofTools"], errors="coerce").fillna(0) >= 2]
    secd = set(l.strip() for l in open(f"{d}/{iso}.secreted_ids.txt") if l.strip())
    def fam(r):
        for c in ("dbCAN_hmm", "DIAMOND", "dbCAN_sub"):
            m_ = re.search(r"(GH\d+|GT\d+|CBM\d+|CE\d+|PL\d+|AA\d+)", str(r[c]))
            if m_:
                return m_.group(1)
    w = ov[ov.apply(fam, axis=1).isin(WALL)]
    return len(w), int(w["Gene ID"].isin(secd).sum())
wc = {i: wall_counts(i) for i in s6.isolate}
s6["Fungal-wall CAZymes (total)"] = s6.isolate.map(lambda i: wc[i][0])
s6["Fungal-wall CAZymes (secreted)"] = s6.isolate.map(lambda i: wc[i][1])
s6["species"] = s6.isolate.map(sp)
s6 = s6.rename(columns={"isolate": "Isolate", "species": "Species",
                        "antifungal_named": "Named antifungal BGC products",
                        "SM:Terpenes (VOCs)": "SM clusters: terpene",
                        "SM:CDPS / pulcherrimin (iron)": "SM clusters: CDPS (pulcherrimin)",
                        "SM:NRPS": "SM clusters: NRPS", "SM:Other": "SM clusters: other"})
s6.columns = [c.replace("b-1,", "β-1,") for c in s6.columns]
s6.to_csv(f"{OUT}/Table_S6_gene_inventory.csv", index=False)

# ---------- Table S7 : killer-toxin tblastn hits ----------
kt = pd.read_csv(f"{D}/killer_toxin_hits.csv")
kt = kt[kt.isolate.isin(SCREENED)]
kt["species"] = kt.isolate.map(sp)
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
# re-deposits of the same strain (see the tsv for ANI evidence) are listed once
dup = set(pd.read_csv("analysis/global_context/meta/duplicate_assemblies.tsv", sep="\t").acc)
s9 = s9[~s9["Assembly accession"].isin(dup)]
s9.to_csv(f"{OUT}/Table_S9_global_comparators.csv", index=False)

# ---------- Table S10 : biocontrol gene content, Kent vs global ----------
bg = pd.read_csv(f"{D}/biocontrol_global_inventory.csv")
# label the public genomes by strain and accession (see analysis/mining/samplesheet_global.tsv)
GLOB = {"JJ10.1": ("JJ10.1", "GCA_026119225.1", "China"),
        "DBVPG 8058": ("DBVPG 8058", "GCA_920103735.2", "unknown (DBVPG collection, Perugia, Italy)"),
        "UKR": ("NRRL Y-1614", "GCA_050230735.1", "Ukraine"),
        "AUS": ("AWRI5759_A6", "GCA_050947715.1", "Australia"),
        "CHN": ("QTX-C10", "GCA_030180015.1", "China")}
bg["accession"] = bg.label.map(lambda l: GLOB[l][1] if l in GLOB else "")
bg["origin"] = [GLOB[l][2] if l in GLOB else o for l, o in zip(bg.label, bg.origin)]
bg["label"] = bg.label.map(lambda l: GLOB[l][0] if l in GLOB else l)
bg["genus"] = bg.genus.replace({"Rhodotorula sp. (aff. babjevae)": "Rhodotorula aff. babjevae"})
bg["grp"] = bg.grp.replace({"global": "Global"})
bg.columns = [c.replace("\n", " ") for c in bg.columns]
cnt = [c for c in bg.columns if c not in ("genus", "label", "origin", "grp", "accession",
                                           "secretome", "n_genes")]
bg = bg[["genus", "label", "accession", "origin", "grp", "n_genes", "secretome"] + cnt]
def _cn(c):
    if c == "SM clusters":
        return "antiSMASH SM clusters"
    return c.split("|")[0] + (" (secreted)" if c.endswith("|sec") else " (total)")
bg.columns = (["Species", "Strain / isolate", "Assembly accession", "Origin", "Group",
               "Predicted genes", "Secreted proteins (SignalP)"] + [_cn(c) for c in cnt])
bg.to_csv(f"{OUT}/Table_S10_biocontrol_kent_vs_global.csv", index=False)

# ---------- summary ----------
for f in sorted(glob.glob(f"{OUT}/Table_S*.csv")):
    df = pd.read_csv(f)
    print(f"{os.path.basename(f):42s} rows={len(df):4d}  cols={len(df.columns)}")
