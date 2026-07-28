#!/usr/bin/env python
"""
31_rename_reads.py — move extracted raw reads into ncbi_upload/reads/ renamed to
paper isolate codes (<isolate>_R1.fastq.gz / _R2.fastq.gz). Records original
Novogene filenames + vendor MD5s in reads_manifest.csv, then removes reads_raw.
"""
import os, glob, shutil, re
import pandas as pd

BASE = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
OUT = f"{BASE}/ncbi_upload"
RAW = f"{OUT}/reads_raw/X204SC24101850-Z01-F001/01.RawData"
READS = f"{OUT}/reads"
os.makedirs(READS, exist_ok=True)

man = pd.read_csv(f"{OUT}/manifest.csv")
sid2iso = dict(zip(man.seq_id, man.isolate))

rows = []
for sid, iso in sid2iso.items():
    d = f"{RAW}/{sid}"
    assert os.path.isdir(d), f"missing extracted dir: {d}"
    # vendor md5s (per-sample MD5.txt: "<md5>  <filename>")
    md5 = {}
    mp = f"{d}/MD5.txt"
    if os.path.exists(mp):
        for ln in open(mp):
            parts = ln.split()
            if len(parts) == 2:
                md5[os.path.basename(parts[1])] = parts[0]
    for rd in ("1", "2"):
        hits = glob.glob(f"{d}/{sid}_*_{rd}.fq.gz")
        assert len(hits) == 1, f"{sid} R{rd}: {hits}"
        src = hits[0]
        orig = os.path.basename(src)
        dst = f"{READS}/{iso}_R{rd}.fastq.gz"
        shutil.move(src, dst)
        rows.append(dict(isolate=iso, seq_id=sid, read=f"R{rd}",
                         new_file=os.path.basename(dst), original_file=orig,
                         vendor_md5=md5.get(orig, "")))

rm = pd.DataFrame(rows).sort_values(["isolate", "read"])
rm.to_csv(f"{OUT}/reads_manifest.csv", index=False)
shutil.rmtree(f"{OUT}/reads_raw")
print(f"renamed {len(rows)} read files for {rm.isolate.nunique()} isolates")
print(f"reads_manifest -> {OUT}/reads_manifest.csv")
