#!/usr/bin/env python
"""
32_clean_assemblies.py — make the renamed assemblies WGS-submission-ready.
For each ncbi_upload/genomes/<isolate>.fasta: drop contigs <200 bp (NCBI auto-trims
these anyway), uppercase, and renumber headers to simple unique IDs (<isolate>_<n>).
Writes cleaned FASTAs to ncbi_upload/genomes_submission/ and a per-genome report.
Originals in genomes/ are left untouched.
"""
import os, glob
import pandas as pd

OUT = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/ncbi_upload"
SRC = f"{OUT}/genomes"
DST = f"{OUT}/genomes_submission"
MINLEN = 200
os.makedirs(DST, exist_ok=True)

def read_fasta(p):
    name, seq = None, []
    for ln in open(p):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = ln[1:], []
        else:
            seq.append(ln)
    if name is not None:
        yield name, "".join(seq)

rows = []
for f in sorted(glob.glob(f"{SRC}/*.fasta")):
    iso = os.path.splitext(os.path.basename(f))[0]
    kept, removed, bp_kept, bp_removed = [], 0, 0, 0
    for _, seq in read_fasta(f):
        s = seq.upper().replace(" ", "")
        if len(s) >= MINLEN:
            kept.append(s)
            bp_kept += len(s)
        else:
            removed += 1
            bp_removed += len(s)
    with open(f"{DST}/{iso}.fasta", "w") as fh:
        for i, s in enumerate(kept, 1):
            fh.write(f">{iso}_{i}\n")
            for j in range(0, len(s), 80):
                fh.write(s[j:j+80] + "\n")
    rows.append(dict(isolate=iso, contigs_kept=len(kept), contigs_removed=removed,
                     bp_kept=bp_kept, bp_removed=bp_removed))

rep = pd.DataFrame(rows).sort_values("isolate")
rep.to_csv(f"{OUT}/assembly_clean_report.csv", index=False)
print(f"cleaned {len(rep)} assemblies -> {DST}")
print(f"  contigs removed (<{MINLEN} bp): total {rep.contigs_removed.sum()}, "
      f"bp removed {rep.bp_removed.sum():,}")
print(f"  report -> {OUT}/assembly_clean_report.csv")
