#!/usr/bin/env python
"""
34_fill_genome_template.py — populate the NCBI batch-genome template (Genome_data sheet)
for the 55 assemblies. Pulls BioSample accessions from the SRA_metadata_acc sheet,
assembler/version from manifest.csv, and genome_coverage from read bases / cleaned
assembly size. Writes Template_GenomeBatch...xlsx in place (openpyxl, formatting kept)
and a parallel .tsv (what actually gets uploaded).
"""
import openpyxl, pandas as pd, glob, os

D = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/ncbi_upload"
GB = f"{D}/Template_GenomeBatch.acdfafce8c89.xlsx"
SRA = f"{D}/SRA_metadata_acc (1).xlsx"

# --- isolate -> BioSample accession (from the registered SRA sheet) ---
sw = openpyxl.load_workbook(SRA, data_only=True)["SRA_data"]
hdr = [c.value for c in sw[1]]
bi, li = hdr.index("biosample_accession"), hdr.index("library_ID")
iso2biosample = {}
for row in sw.iter_rows(min_row=2, values_only=True):
    if row[li]:
        iso2biosample[str(row[li])] = row[bi]
assert len(iso2biosample) == 55, len(iso2biosample)

# --- assembler/version + assembly date ---
man = pd.read_csv(f"{D}/manifest.csv").set_index("isolate")
def method(iso):
    a = man.loc[iso, "assembler"]              # e.g. "SPAdes v3.15.5"
    if a.startswith("Platanus"):
        return "Platanus-allee", "2.4.0", "2026-07"
    return "SPAdes", "3.15.5", "2024-12"

# --- genome coverage (read bases / cleaned assembly size) ---
# Built by scratchpad/coverage.py: 41 isolates from exact base counts, 14 predicted
# from compressed read size (bases/byte CV 3.1%); validated vs mapping_qc depth to 1.1%.
cov = pd.read_csv("/tmp/coverage.csv").set_index("isolate")["coverage"]
assert set(cov.index) == set(iso2biosample), \
    f"isolate mismatch: coverage={len(cov)} biosamples={len(iso2biosample)}"

wb = openpyxl.load_workbook(GB)
ws = wb["Genome_data"]
cols = [c.value for c in ws[1]]   # header order
def put(r, name, val):
    ws.cell(row=r, column=cols.index(name) + 1, value=val)

rows = []
for i, iso in enumerate(sorted(iso2biosample), start=2):
    meth, ver, adate = method(iso)
    rec = dict(biosample_accession=iso2biosample[iso], sample_name=iso,
               assembly_date=adate, assembly_name=f"{iso}_v1.0",
               assembly_method=meth, assembly_method_version=ver,
               genome_coverage=f"{int(cov.loc[iso])}x", sequencing_technology="Illumina NovaSeq X Plus",
               reference_genome="", update_for="", bacteria_available_from="",
               filename=f"{iso}.fasta")
    for k, v in rec.items():
        put(i, k, v)
    rows.append(rec)

wb.save(GB)
hdr_cols = [c for c in cols if c]           # drop trailing empty header cells
pd.DataFrame(rows)[hdr_cols].to_csv(f"{D}/genome_data.tsv", sep="\t", index=False)
print(f"filled {len(rows)} rows -> {os.path.basename(GB)} (Genome_data) + genome_data.tsv")
cv = [int(r["genome_coverage"][:-1]) for r in rows]
print(f"coverage range: {min(cv)}x – {max(cv)}x  (median {sorted(cv)[len(cv)//2]}x)")
print(f"methods: SPAdes {sum(r['assembly_method']=='SPAdes' for r in rows)}, "
      f"Platanus-allee {sum(r['assembly_method']=='Platanus-allee' for r in rows)}")
