#!/usr/bin/env python
"""
33_ncbi_metadata.py — pre-fill NCBI submission spreadsheets from the isolate master.
Writes biosample_attributes.tsv (Microbe 1.0 package) and sra_metadata.tsv for the
55 paper isolates. 'organism' uses an NCBI-valid name (genus 'sp.' where our genome
ID is a 'cf.'/undescribed lineage); the precise determination is kept in identified_as.
"""
import pandas as pd

BASE = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
OUT = f"{BASE}/ncbi_upload"

man = pd.read_csv(f"{OUT}/manifest.csv")   # 55, has isolate/seq_id/species

# genome determination -> NCBI-submittable organism name (cf./undescribed -> genus sp.)
ORG = {
    "Hanseniaspora uvarum": "Hanseniaspora uvarum",
    "Pichia kluyveri": "Pichia kluyveri",
    "Sporobolomyces roseus": "Sporobolomyces roseus",
    "Sporobolomyces sp.": "Sporobolomyces sp.",
    "Rhodotorula aff. babjevae": "Rhodotorula sp.",
    "Metschnikowia cf. pulcherrima": "Metschnikowia sp.",
    "Vishniacozyma carnescens": "Vishniacozyma carnescens",
    "Vishniacozyma victoriae": "Vishniacozyma victoriae",
    "Vishniacozyma cf. victoriae": "Vishniacozyma sp.",
    "Filobasidium sp. (cf. globisporum)": "Filobasidium sp.",
}

def source(iso):
    if iso.startswith("UY"):
        return ("Vitis vinifera", "Chardonnay",
                "surface of grape berry (Vitis vinifera cv. Chardonnay), untreated vineyard block")
    if iso.startswith("TY"):
        return ("Vitis vinifera", "Chardonnay",
                "surface of grape berry (Vitis vinifera cv. Chardonnay), fungicide-treated vineyard block")
    return ("Vitis vinifera", "Pinot Noir",
            "surface of grape berry (Vitis vinifera cv. Pinot Noir)")

GEO = "United Kingdom: England, Kent, East Malling"
LATLON = "51.29 N 0.45 E"      # NIAB East Malling research vineyard (approximate)
CDATE = "2024-09"

bio, sra = [], []
for _, r in man.sort_values("isolate").iterrows():
    iso = r.isolate
    host, cv, isrc = source(iso)
    bio.append({
        "sample_name": iso,
        "organism": ORG[r.species],
        "isolate": iso,
        "isolation_source": isrc,
        "collection_date": CDATE,
        "geo_loc_name": GEO,
        "lat_lon": LATLON,
        "host": host,
        "host_cultivar": cv,
        "identified_as": r.species,            # full genome-based determination
    })
    sra.append({
        "sample_name": iso,
        "library_ID": iso,
        "title": f"Whole-genome sequencing of vineyard yeast isolate {iso} ({r.species})",
        "library_strategy": "WGS",
        "library_source": "GENOMIC",
        "library_selection": "RANDOM",
        "library_layout": "paired",
        "platform": "ILLUMINA",
        "instrument_model": "Illumina NovaSeq X Plus",
        "design_description": "PCR-free/short-insert paired-end library, PE150 (Novogene)",
        "filetype": "fastq",
        "filename": f"{iso}_R1.fastq.gz",
        "filename2": f"{iso}_R2.fastq.gz",
    })

pd.DataFrame(bio).to_csv(f"{OUT}/biosample_attributes.tsv", sep="\t", index=False)
pd.DataFrame(sra).to_csv(f"{OUT}/sra_metadata.tsv", sep="\t", index=False)
print(f"biosample_attributes.tsv  ({len(bio)} rows)")
print(f"sra_metadata.tsv          ({len(sra)} rows)")
print("\norganism counts:")
print(pd.DataFrame(bio).organism.value_counts().to_string())
