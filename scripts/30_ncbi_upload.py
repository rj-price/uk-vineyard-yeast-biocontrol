#!/usr/bin/env python
"""
30_ncbi_upload.py — assemble the NCBI submission package.
Copies genome assemblies (renamed to paper isolate codes) into ncbi_upload/genomes/,
writes a manifest, and emits the tar member-list for selective read extraction.
Paper set = 55 sequenced isolates (excludes 5 dropped Sporobolomyces).
"""
import os, glob, shutil
import pandas as pd

BASE = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
OUT = f"{BASE}/ncbi_upload"
GEN = f"{OUT}/genomes"
PLAT = f"{BASE}/analysis/final_assemblies/ascomycetes_platanus"   # 26 ascomycete Platanus
SPADES = f"{BASE}/genomics"                                        # SPAdes scaffolds
TARROOT = "X204SC24101850-Z01-F001/01.RawData"

DROP = {"TY8", "TY9", "TY12", "TY17", "TY24"}   # 5 dropped Sporobolomyces
os.makedirs(GEN, exist_ok=True)

m = pd.read_csv(f"{BASE}/analysis/data/isolate_master.csv")
m = m[~m.isolate.isin(DROP)].copy()
assert len(m) == 55, f"expected 55, got {len(m)}"

# ascomycete = has a Platanus assembly named by seq_id
plat_ids = {os.path.splitext(os.path.basename(p))[0] for p in glob.glob(f"{PLAT}/*.fasta")}

rows, members = [], []
for _, r in m.iterrows():
    iso, sid, sp = r.isolate, r.seq_id, r.species_genome
    if sid in plat_ids:
        src = f"{PLAT}/{sid}.fasta"
        asm = "Platanus-allee v2.4.0"
    else:
        hits = glob.glob(f"{SPADES}/{sid}_*_scaffolds.fasta")
        assert len(hits) == 1, f"{sid}: {hits}"
        src = hits[0]
        asm = "SPAdes v3.15.5"
    dst = f"{GEN}/{iso}.fasta"
    shutil.copyfile(src, dst)
    # read files: barcode string is embedded in the SPAdes filename; derive it,
    # else fall back to extracting the whole seq_id dir from the tar.
    members.append(f"{TARROOT}/{sid}")
    rows.append(dict(isolate=iso, seq_id=sid, species=sp, assembler=asm,
                     assembly_source=os.path.relpath(src, BASE)))

man = pd.DataFrame(rows).sort_values("isolate")
man.to_csv(f"{OUT}/manifest.csv", index=False)
with open(f"{OUT}/tar_members.txt", "w") as fh:
    fh.write("\n".join(members) + "\n")

print(f"genomes copied: {len(rows)}")
print(f"  Platanus: {(man.assembler.str.startswith('Platanus')).sum()}  "
      f"SPAdes: {(man.assembler.str.startswith('SPAdes')).sum()}")
print(f"manifest -> {OUT}/manifest.csv")
print(f"members  -> {OUT}/tar_members.txt ({len(members)} dirs)")
