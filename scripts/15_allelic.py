#!/usr/bin/env python
"""
15_allelic.py
Allelic + upstream variation in Hanseniaspora uvarum, across the B. cinerea
antagonism gradient RY14 (62%) > RY7 (36.5%) > RY25 (13.7%).

NOTE ON FRAMING: the dual-plate replicate data show all three are equally
inactive against BOTH Fusaria (85/85/85 = plate ceiling, every replicate), so
the only trait that varies here is B. cinerea suppression. Gene-content mining
(11_comparative.py) found no presence/absence difference, so we ask instead
whether biocontrol genes carry RY14-specific coding or upstream (putative
promoter) variants at a rate above the genome-wide background.

CRITICAL: RY14 is a single genome, so "RY14-specific" variants accumulate
wherever it is simply more diverged. Within-group ANI is flat (98.85-99.0%,
no structure), so RY14 is not an outlier lineage - but the genome-wide
background distribution remains the only meaningful reference: a biocontrol
gene is interesting only if it is an OUTLIER against that background.

Rhodotorula is handled separately (03_clonal_pair.sh): UY1/TY27 are 99.99% ANI
clonal, so their variants are enumerated exhaustively rather than tested.

-> analysis/data/allelic_<genus>.csv        (per-orthogroup stats)
-> analysis/data/allelic_summary.csv        (biocontrol vs background test)
"""
import os, re, subprocess, sys, tempfile
import pandas as pd, numpy as np
from Bio import SeqIO
from scipy import stats

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
ALL = f"{PROJ}/analysis/allelic"
RES = f"{PROJ}/analysis/mining/results"
DATA = f"{PROJ}/analysis/data"
# call the binary directly: `conda run` costs ~1.5s of startup per invocation,
# which over thousands of orthogroups x 2 alignments would dominate the runtime
ENV = "/mnt/apps/users/jnprice/conda/envs/orthofinder"
MAFFT = f"{ENV}/bin/mafft"
MAFFT_ENV = dict(os.environ, MAFFT_BINARIES=f"{ENV}/libexec/mafft")
UPSTREAM = 1000

GROUPS = {
    "Hanseniaspora": dict(antag="RY14", others=["RY7", "RY25"]),
}


def load_annotations(iso):
    """Return (transcript -> annotation row) and the secreted transcript set."""
    a = pd.read_csv(f"{RES}/{iso}/{iso}.annotations.txt", sep="\t", dtype=str).fillna("")
    sec = set(open(f"{RES}/{iso}/{iso}.secreted_ids.txt").read().split())
    return a, sec


def biocontrol_transcripts(iso):
    """Transcripts of the antagonist that are secreted AND CAZyme/protease-annotated."""
    a, sec = load_annotations(iso)
    hit = a[(a.TranscriptID.isin(sec)) & ((a.CAZyme != "") | (a.Protease != ""))]
    return dict(zip(hit.TranscriptID, hit.CAZyme + "|" + hit.Protease + "|" + hit.Product))


def read_orthogroups(genus):
    """Single-copy orthogroups: exactly one transcript per isolate."""
    base = f"{ALL}/{genus}/OrthoFinder"
    run = sorted(os.listdir(base))[-1]
    og = pd.read_csv(f"{base}/{run}/Orthogroups/Orthogroups.tsv", sep="\t", dtype=str)
    og = og.set_index("Orthogroup")
    single = og[og.notna().all(axis=1) & og.apply(
        lambda r: all("," not in str(v) for v in r), axis=1)]
    return single


def load_prots(genus, isos):
    d = {}
    for i in isos:
        d[i] = {r.id: str(r.seq).rstrip("*") for r in
                SeqIO.parse(f"{ALL}/{genus}/{i}.fa", "fasta")}
    return d


def gene_coords(iso):
    """gene ID -> (contig, start, stop, strand) from the funannotate gff3."""
    out = {}
    for line in open(f"{RES}/{iso}/{iso}.gff3"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        gid = re.search(r"ID=([^;]+)", f[8]).group(1)
        out[gid] = (f[0], int(f[3]), int(f[4]), f[6])
    return out


def load_assembly(iso, subjects=None):
    """MUST use funannotate's sorted.fa, not the original assembly: funannotate
    renames contigs (NODE_*/platanus -> scaffold_*) during sort, so the gff3
    coordinates only resolve against sorted.fa."""
    path = f"{PROJ}/analysis/mining/work/{iso}/sorted.fa"
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def upstream_seq(gid, coords, contigs):
    if gid not in coords:
        return None
    ctg, start, stop, strand = coords[gid]
    if ctg not in contigs:
        return None
    s = contigs[ctg]
    if strand == "+":
        seq = s[max(0, start - 1 - UPSTREAM): start - 1]
    else:
        seq = s[stop: stop + UPSTREAM]
        seq = seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    return seq if len(seq) >= 100 else None


def mafft_align(seqs, nucl=False):
    """seqs: {name: seq} -> {name: aligned}. Returns None on failure."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n{s}\n")
        tmp = fh.name
    cmd = [MAFFT, "--auto", "--quiet", "--nuc" if nucl else "--amino", tmp]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                           env=MAFFT_ENV)
        if p.returncode != 0 or not p.stdout.strip():
            return None
        out, name = {}, None
        for line in p.stdout.splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = ""
            elif name:
                out[name] += line.strip()
        return {k: v.upper() for k, v in out.items()}
    except subprocess.TimeoutExpired:
        return None
    finally:
        os.unlink(tmp)


def specific_sites(aln, antag, others):
    """Count alignment columns where the antagonist differs from ALL others and
    the others agree with each other (a clean antagonist-specific state), plus
    total variable columns. Gapped columns are skipped."""
    if aln is None or antag not in aln or any(o not in aln for o in others):
        return None
    L = len(aln[antag])
    if any(len(aln[o]) != L for o in others):
        return None
    spec = var = usable = 0
    for k in range(L):
        a = aln[antag][k]
        os_ = [aln[o][k] for o in others]
        if a == "-" or any(x == "-" for x in os_):
            continue
        usable += 1
        if len(set(os_ + [a])) > 1:
            var += 1
        if len(set(os_)) == 1 and a != os_[0]:
            spec += 1
    return dict(len_aln=usable, n_var=var, n_spec=spec)


def run_genus(genus, subjects):
    antag = GROUPS[genus]["antag"]
    others = GROUPS[genus]["others"]
    isos = [antag] + others
    print(f"\n### {genus}: antagonist {antag} vs {others}", flush=True)

    single = read_orthogroups(genus)
    print(f"  single-copy orthogroups: {len(single)}", flush=True)
    prots = load_prots(genus, isos)
    bc = biocontrol_transcripts(antag)
    coords = {i: gene_coords(i) for i in isos}
    contigs = {i: load_assembly(i) for i in isos}
    for i in isos:  # fail loudly rather than silently emitting all-NaN upstream
        hit = sum(1 for c, *_ in coords[i].values() if c in contigs[i])
        if hit == 0:
            raise SystemExit(f"{i}: no gff3 contig matches sorted.fa - naming mismatch")
        print(f"  {i}: {hit}/{len(coords[i])} genes on loaded contigs", flush=True)

    rows = []
    for n, (og, r) in enumerate(single.iterrows()):
        tx = {i: str(r[i]).split()[0] for i in isos}
        pseqs = {i: prots[i].get(tx[i]) for i in isos}
        if any(s is None for s in pseqs.values()):
            continue
        st = specific_sites(mafft_align(pseqs), antag, others)
        if st is None or st["len_aln"] < 50:
            continue

        useqs = {i: upstream_seq(tx[i].rsplit("-T", 1)[0], coords[i], contigs[i])
                 for i in isos}
        ust = (specific_sites(mafft_align(useqs, nucl=True), antag, others)
               if all(u for u in useqs.values()) else None)

        rows.append(dict(
            orthogroup=og, antag_tx=tx[antag],
            aa_len=st["len_aln"], aa_var=st["n_var"], aa_spec=st["n_spec"],
            up_len=ust["len_aln"] if ust else np.nan,
            up_var=ust["n_var"] if ust else np.nan,
            up_spec=ust["n_spec"] if ust else np.nan,
            biocontrol=tx[antag] in bc, bc_annot=bc.get(tx[antag], "")))
        if (n + 1) % 500 == 0:
            print(f"    ...{n+1} orthogroups processed", flush=True)

    d = pd.DataFrame(rows)
    d["aa_spec_rate"] = d.aa_spec / d.aa_len
    d["up_spec_rate"] = d.up_spec / d.up_len
    d.to_csv(f"{DATA}/allelic_{genus}.csv", index=False)
    return d


def summarise(genus, d):
    out = []
    for lab, rate in [("protein", "aa_spec_rate"), ("upstream", "up_spec_rate")]:
        sub = d.dropna(subset=[rate])
        b = sub[sub.biocontrol][rate]
        g = sub[~sub.biocontrol][rate]
        if len(b) < 3:
            continue
        u = stats.mannwhitneyu(b, g, alternative="two-sided")
        out.append(dict(genus=genus, level=lab, n_biocontrol=len(b), n_background=len(g),
                        median_bc=round(b.median(), 5), median_bg=round(g.median(), 5),
                        p_mannwhitney=round(u.pvalue, 4)))
    return out


def main():
    os.chdir(PROJ)
    subjects = pd.read_csv("analysis/mining/killer_subjects.tsv", sep="\t")
    summ = []
    for genus in GROUPS:
        d = run_genus(genus, subjects)
        summ += summarise(genus, d)

        bg = d[~d.biocontrol].aa_spec_rate
        thr = bg.quantile(0.95)
        print(f"\n  genome-wide background antagonist-specific aa rate: "
              f"median {bg.median():.4f}, 95th pct {thr:.4f}")
        outl = d[d.biocontrol & (d.aa_spec_rate > thr)].sort_values(
            "aa_spec_rate", ascending=False)
        print(f"  biocontrol genes above the 95th-pct background "
              f"({len(outl)}/{d.biocontrol.sum()}):")
        if len(outl):
            print(outl[["antag_tx", "aa_len", "aa_spec", "aa_spec_rate",
                        "bc_annot"]].head(15).to_string(index=False))

    s = pd.DataFrame(summ)
    s.to_csv(f"{DATA}/allelic_summary.csv", index=False)
    print("\n=== biocontrol vs genome-wide background ===")
    print(s.to_string(index=False))
    print(f"\nWrote {DATA}/allelic_<genus>.csv + allelic_summary.csv")


if __name__ == "__main__":
    main()
