#!/usr/bin/env python
"""
16_clonal_variants.py
Annotate the complete variant list separating the clonal pair UY1 (B. cinerea+,
34.5% MI) and TY27 (B. cinerea-, 0% MI), 99.99% ANI. Variants are the read-validated
calls from 38_clonal_reads.py (the earlier assembly-vs-assembly dnadiff set was ~86%
assembly artefact).

Because divergence is so low, no statistical outlier test is needed or possible:
the variant list IS the candidate set. We classify every variant by gene context
and effect, and flag those touching biocontrol-relevant genes.

TY28/TY29 (~98.3% ANI) are NOT useful outgroups at SNP level (~265k variants each,
so they differ nearly everywhere); they are used only to check whether the UY1
allele at a candidate site is the derived/rare one.

CAVEAT: n=1 clonal pair. Any candidate is a correlation, not a demonstrated cause.

-> analysis/data/clonal_variants.csv        (every variant, annotated)
-> analysis/data/clonal_candidates.csv      (variants in/near biocontrol genes)
"""
import os, re
import pandas as pd, numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
CL = f"{PROJ}/analysis/allelic/clonal_ref"
RES = f"{PROJ}/analysis/mining/results"
DATA = f"{PROJ}/analysis/data"
REF, QRY = "UY1", "TY27"
UPSTREAM = 1000


def read_snps(path):
    """dnadiff .snps (show-snps -rlTH): rpos rbase qbase qpos ... rcontig qcontig.
    A '.' in rbase/qbase marks an indel."""
    rows = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        rows.append(dict(rpos=int(f[0]), rbase=f[1], qbase=f[2], qpos=int(f[3]),
                         rctg=f[-2], qctg=f[-1]))
    d = pd.DataFrame(rows)
    d["kind"] = np.where((d.rbase == ".") | (d.qbase == "."), "indel", "snp")
    return d


def load_genes(iso):
    """mRNA features with CDS parts, from the funannotate gff3."""
    genes, cds = {}, {}
    for line in open(f"{RES}/{iso}/{iso}.gff3"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        att = f[8]
        if f[2] == "mRNA":
            tid = re.search(r"ID=([^;]+)", att).group(1)
            genes[tid] = dict(ctg=f[0], start=int(f[3]), stop=int(f[4]), strand=f[6])
        elif f[2] == "CDS":
            pid = re.search(r"Parent=([^;]+)", att).group(1)
            cds.setdefault(pid, []).append((int(f[3]), int(f[4]), int(f[7]) if f[7].isdigit() else 0))
    for t in cds:
        cds[t].sort()
    return genes, cds


def annot_table(iso):
    a = pd.read_csv(f"{RES}/{iso}/{iso}.annotations.txt", sep="\t", dtype=str).fillna("")
    sec = set(open(f"{RES}/{iso}/{iso}.secreted_ids.txt").read().split())
    a["is_secreted"] = a.TranscriptID.isin(sec)
    return a.set_index("TranscriptID")


def build_index(genes):
    """contig -> list of (start, stop, tid) for genes and their upstream windows."""
    gi, ui = {}, {}
    for t, g in genes.items():
        gi.setdefault(g["ctg"], []).append((g["start"], g["stop"], t))
        if g["strand"] == "+":
            u = (max(1, g["start"] - UPSTREAM), g["start"] - 1, t)
        else:
            u = (g["stop"] + 1, g["stop"] + UPSTREAM, t)
        ui.setdefault(g["ctg"], []).append(u)
    return gi, ui


def locate(ctg, pos, idx):
    return [t for s, e, t in idx.get(ctg, []) if s <= pos <= e]


def cds_effect(tid, pos, rbase, qbase, genes, cds, contigs):
    """Return (effect, aa_change, codon_pos) for a SNP inside a CDS."""
    if tid not in cds:
        return ("noncoding_exon", "", np.nan)
    g = genes[tid]
    parts = cds[tid]
    seq = "".join(contigs[g["ctg"]][s - 1:e] for s, e, _ in parts)
    # offset of pos within the concatenated CDS (genomic order)
    off, found = 0, None
    for s, e, _ in parts:
        if s <= pos <= e:
            found = off + (pos - s)
            break
        off += e - s + 1
    if found is None:
        return ("intron", "", np.nan)
    ref = Seq(seq)
    alt = Seq(seq[:found] + qbase + seq[found + 1:])
    if g["strand"] == "-":
        ref, alt = ref.reverse_complement(), alt.reverse_complement()
        found = len(seq) - 1 - found
    ci = found // 3
    try:
        pr, pa = str(ref.translate()), str(alt.translate())
    except Exception:
        return ("cds_untranslatable", "", np.nan)
    if ci >= len(pr) or ci >= len(pa):
        return ("cds_out_of_frame", "", np.nan)
    if pr[ci] == pa[ci]:
        return ("synonymous", "", ci + 1)
    eff = "nonsense" if pa[ci] == "*" else ("stop_lost" if pr[ci] == "*" else "missense")
    return (eff, f"{pr[ci]}{ci+1}{pa[ci]}", ci + 1)


def read_validated():
    """Read-validated UY1/TY27 differences from 38_clonal_reads.py (VCF-style alleles)."""
    v = pd.read_csv(f"{DATA}/clonal_reads_validated.csv")
    v = v[(v.group == REF) & (v[f"{QRY}_gt"] > 0)]
    return pd.DataFrame(dict(rpos=v.pos, rbase=v.ref, qbase=v.alt.str.split(",").str[0],
                             rctg=v.ctg, kind=np.where(v.indel, "indel", "snp")))


def indel_effect(tid, pos, rbase, qbase, genes, cds):
    """Frameshift / in-frame for an indel overlapping a CDS, else intron/UTR."""
    if tid not in cds:
        return "indel_noncoding_exon"
    span = range(pos, pos + len(rbase))
    if not any(s <= p <= e for s, e, _ in cds[tid] for p in span):
        return "indel_intron_or_utr"
    return "inframe_indel" if (len(qbase) - len(rbase)) % 3 == 0 else "frameshift"


def main():
    os.chdir(PROJ)
    # read-level calls (reviewer point 3); the dnadiff set is kept for comparison in
    # clonal_variants_dnadiff.csv / clonal_dnadiff_vs_reads.csv
    d = read_validated()
    print(f"variants UY1 vs {QRY}: {len(d)} "
          f"({(d.kind=='snp').sum()} SNPs, {(d.kind=='indel').sum()} indel rows)")

    genes, cds = load_genes(REF)
    contigs = {r.id: str(r.seq).upper()
               for r in SeqIO.parse(f"{PROJ}/analysis/mining/work/{REF}/sorted.fa", "fasta")}
    ann = annot_table(REF)
    gi, ui = build_index(genes)

    rows = []
    for r in d.itertuples():
        ingene = locate(r.rctg, r.rpos, gi)
        inup = locate(r.rctg, r.rpos, ui)
        if ingene:
            tid = ingene[0]
            ctx = "genic"
            eff, aac, cpos = (cds_effect(tid, r.rpos, r.rbase, r.qbase, genes, cds, contigs)
                             if r.kind == "snp" else
                             (indel_effect(tid, r.rpos, r.rbase, r.qbase, genes, cds), "", np.nan))
        elif inup:
            tid, ctx, eff, aac, cpos = inup[0], "upstream", "upstream", "", np.nan
        else:
            tid, ctx, eff, aac, cpos = "", "intergenic", "intergenic", "", np.nan
        a = ann.loc[tid] if tid in ann.index else None
        rows.append(dict(
            contig=r.rctg, pos=r.rpos, ref=r.rbase, alt=r.qbase, kind=r.kind,
            context=ctx, effect=eff, aa_change=aac, codon=cpos, transcript=tid,
            product=(a.Product if a is not None else ""),
            cazyme=(a.CAZyme if a is not None else ""),
            protease=(a.Protease if a is not None else ""),
            secreted=(bool(a.is_secreted) if a is not None else False),
            pfam=(a.PFAM if a is not None else "")))
    v = pd.DataFrame(rows)
    v.to_csv(f"{DATA}/clonal_variants.csv", index=False)

    print("\n=== variants by context ===")
    print(v.context.value_counts().to_string())
    print("\n=== coding effects (SNPs in genes) ===")
    print(v[v.context == "genic"].effect.value_counts().to_string())

    v["biocontrol"] = (v.cazyme != "") | (v.protease != "") | v.secreted
    cand = v[v.biocontrol & v.context.isin(["genic", "upstream"])].copy()
    cand.to_csv(f"{DATA}/clonal_candidates.csv", index=False)

    print(f"\n=== variants touching biocontrol genes (CAZyme/protease/secreted): "
          f"{len(cand)} in {cand.transcript.nunique()} genes ===")
    if len(cand):
        show = cand[["transcript", "context", "effect", "aa_change", "cazyme",
                     "protease", "secreted", "product"]]
        print(show.to_string(index=False))

    imp = v[v.effect.isin(["missense", "nonsense", "stop_lost", "frameshift", "inframe_indel"])]
    print(f"\n=== ALL protein-altering variants genome-wide: {len(imp)} in "
          f"{imp.transcript.nunique()} genes ===")
    print(imp.effect.value_counts().to_string())
    print(f"\nWrote {DATA}/clonal_variants.csv + clonal_candidates.csv")


if __name__ == "__main__":
    main()
