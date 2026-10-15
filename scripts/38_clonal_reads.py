#!/usr/bin/env python
"""
38_clonal_reads.py
Read-level validation of within-clone variants (reviewer point 3), from the joint
haploid calls of allelic/07_clonal_call.sh.

Filters (per site): QUAL >= 30; every member DP >= 10 and called; the called allele
carries >= 90% of that member's reads (haploid; mixed sites usually mean collapsed
paralogues). Sites where the reference isolate's own reads call ALT are assembly
errors in the reference and are counted separately.

For group A, the assembly-based dnadiff set (data/clonal_variants_dnadiff.csv, UY1 vs
TY27 on the same UY1 scaffolds) is intersected with the read-validated UY1/TY27 differences.

-> analysis/data/clonal_reads_summary.csv
-> analysis/data/clonal_reads_validated.csv   (filtered between-member variants; input to 16)
-> analysis/data/clonal_dnadiff_vs_reads.csv  (dnadiff variants flagged read_confirmed)
"""
import gzip, os
import numpy as np, pandas as pd

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis"
CR = f"{PROJ}/allelic/clonal_reads"
GROUPS = {"UY1": ["UY1", "TY27", "UY4"], "TY28": ["TY28", "UY2", "UY3"]}


def read_vcf(path):
    rows = []
    for line in gzip.open(path, "rt"):
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            samples = f[9:]; continue
        fmt = f[8].split(":")
        rec = dict(ctg=f[0], pos=int(f[1]), ref=f[3], alt=f[4], qual=float(f[5]) if f[5] != "." else 0.0,
                   indel=len(f[3]) != 1 or any(len(a) != 1 for a in f[4].split(",")))
        for s, v in zip(samples, f[9:]):
            d = dict(zip(fmt, v.split(":")))
            gt = d.get("GT", ".")
            ad = [int(x) for x in d.get("AD", "0").split(",") if x != "."]
            dp = sum(ad)
            g = int(gt) if gt not in (".", "./.") else -1
            rec[s + "_gt"] = g
            rec[s + "_dp"] = dp
            rec[s + "_af"] = (ad[g] / dp) if (g >= 0 and dp and g < len(ad)) else 0.0
        rows.append(rec)
    return pd.DataFrame(rows)


rows, keep = [], []
for ref, mem in GROUPS.items():
    if not os.path.getsize(f"{CR}/group_{ref}.vcf.gz"):
        print(f"group {ref}: VCF not ready, skipped"); continue
    v = read_vcf(f"{CR}/group_{ref}.vcf.gz")
    ok = (v.qual >= 30)
    for s in mem:
        ok &= (v[s + "_gt"] >= 0) & (v[s + "_dp"] >= 10) & (v[s + "_af"] >= 0.9)
    f = v[ok].copy()
    self_err = f[f[ref + "_gt"] > 0]
    f = f[f[ref + "_gt"] == 0]                         # reference reads agree with the reference
    for s in mem[1:]:
        d = f[f[s + "_gt"] > 0]
        rows.append(dict(group=ref, pair=f"{ref} vs {s}", raw_sites=len(v), pass_filters=int(ok.sum()),
                         ref_self_errors=len(self_err), variants=len(d),
                         snps=int((~d.indel).sum()), indels=int(d.indel.sum())))
    a, b = mem[1], mem[2]
    sh = f[(f[a + "_gt"] > 0) & (f[b + "_gt"] > 0)]
    rows.append(dict(group=ref, pair=f"shared by {a} and {b} (private to {ref})", variants=len(sh),
                     snps=int((~sh.indel).sum()), indels=int(sh.indel.sum())))
    f.insert(0, "group", ref)
    keep.append(f)

summ = pd.DataFrame(rows)
val = pd.concat(keep, ignore_index=True)

# dnadiff (assembly) UY1 vs TY27 set, checked against reads
dn = pd.read_csv(f"{PROJ}/data/clonal_variants_dnadiff.csv")   # frozen dnadiff-based table
a = val[(val.group == "UY1") & (val.TY27_gt > 0)]
pos = set(zip(a.ctg, a.pos))
# indels can be left-shifted differently: match within 5 bp on the same scaffold
by_ctg = a.groupby("ctg").pos.apply(np.array).to_dict()
hit = [any(abs(by_ctg.get(c, np.array([-99]))[:, None] - p).min() <= 5 for _ in [0]) for c, p in zip(dn.contig, dn.pos)]
dn["read_confirmed"] = hit
for k, g in dn.groupby("kind"):
    print(f"dnadiff UY1/TY27 {k}: {len(g)} sites, {int(g.read_confirmed.sum())} confirmed by reads")
pa = dn[dn.effect.isin(["missense", "nonsense", "indel_in_gene", "frameshift", "stop_lost", "start_lost"])]
print(f"  protein-altering (Table S8 set): {len(pa)}, confirmed {int(pa.read_confirmed.sum())}")
for t in ["UY1_007378-T1", "UY1_003568-T1"]:
    g = dn[dn.transcript == t]
    print(f"  {t}: {len(g)} dnadiff variants, {int(g.read_confirmed.sum())} confirmed")
summ = pd.concat([summ, pd.DataFrame([dict(group="UY1", pair="dnadiff UY1 vs TY27 confirmed by reads",
                                           variants=int(dn.read_confirmed.sum()), raw_sites=len(dn))])])
summ.to_csv(f"{PROJ}/data/clonal_reads_summary.csv", index=False)
val.to_csv(f"{PROJ}/data/clonal_reads_validated.csv", index=False)
dn.to_csv(f"{PROJ}/data/clonal_dnadiff_vs_reads.csv", index=False)
print(summ.to_string(index=False))
