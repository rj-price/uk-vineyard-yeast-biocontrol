#!/usr/bin/env python
"""
18_interpro_enrichment.py
Repeat the missense/protein-altering enrichment for the UY1(Bc+)/TY27(Bc-) clonal pair
using proper InterPro + GO annotation instead of PFAM-regex.

The statistics are NOT expected to change: the binding constraint is the gene-length
confound (regulators are long; long genes accrue more variants), not annotation depth.
The purpose is to surface better-characterised candidate loci for Discussion, and to
replace the earlier hand-written PFAM regex (which had false positives, e.g. it classed
a trehalose-6-phosphate phosphatase as "signalling") with curated GO terms.

-> analysis/data/interpro_annotation.csv    (per-gene IPR/GO, all 12 genomes)
-> analysis/data/clonal_variants_ipr.csv    (the 97 variant-bearing genes, annotated)
-> analysis/data/interpro_enrichment.csv    (GO-category enrichment, length-controlled)
"""
import os, glob
import pandas as pd, numpy as np
from scipy import stats
from Bio import SeqIO
import statsmodels.api as sm

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
IPRO = f"{PROJ}/analysis/interpro/results"
DATA = f"{PROJ}/analysis/data"
REF = os.environ.get("CLONAL_REF", "UY1")
SFX = "" if REF == "UY1" else f"_{REF}"

# Curated GO terms per category (ancestor terms are not propagated by InterProScan,
# so list the specific terms that actually appear on domain matches).
GO_CATS = {
    "transcription_regulation": [
        "GO:0006355",  # regulation of DNA-templated transcription
        "GO:0003700",  # DNA-binding transcription factor activity
        "GO:0000981",  # DNA-binding TF activity, RNA pol II-specific
        "GO:0043565",  # sequence-specific DNA binding
        "GO:0000122", "GO:0045893", "GO:0045892",
        "GO:0006357",  # regulation of transcription by RNA pol II
    ],
    "chromatin_epigenetic": [
        "GO:0006325",  # chromatin organization
        "GO:0016570",  # histone modification
        "GO:0018024", "GO:0042054",  # histone methyltransferase activity
        "GO:0004402",  # histone acetyltransferase activity
        "GO:0031507",  # heterochromatin formation
        "GO:0006306", "GO:0003682",  # DNA methylation, chromatin binding
        "GO:0032259",  # methylation
    ],
    "signal_transduction": [
        "GO:0007165",  # signal transduction
        "GO:0035556",  # intracellular signal transduction
        "GO:0004672",  # protein kinase activity
        "GO:0006468",  # protein phosphorylation
        "GO:0004674", "GO:0004713",  # Ser/Thr and Tyr kinase activity
        "GO:0007124", "GO:0000165",  # pseudohyphal growth, MAPK cascade
        # NB: deliberately NOT GO:0005515 ("protein binding") - it is the most generic
        # term in GO and hit 507/4641 genes in testing, i.e. pure noise, not signalling.
    ],
}


def read_ips(path):
    cols = ["protein", "md5", "length", "analysis", "sig_acc", "sig_desc",
            "start", "stop", "evalue", "status", "date", "ipr_acc", "ipr_desc", "go", "pathway"]
    d = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str,
                    on_bad_lines="skip").fillna("")
    return d


def per_gene(d):
    """Collapse InterProScan hits to one row per protein."""
    d = d[d.ipr_acc.str.startswith("IPR")]
    g = d.groupby("protein").agg(
        ipr=("ipr_acc", lambda s: "|".join(sorted(set(s)))),
        ipr_desc=("ipr_desc", lambda s: "; ".join(sorted(set(s)))),
        go=("go", lambda s: "|".join(sorted({x.split("(")[0] for v in s for x in v.split("|")
                                             if x.startswith("GO:")}))),
    ).reset_index()
    return g


def main():
    os.chdir(PROJ)
    files = sorted(glob.glob(f"{IPRO}/*.tsv"))
    if not files:
        raise SystemExit(f"no InterProScan results in {IPRO} - run 01_interproscan.sh first")
    print(f"InterProScan results found: {len(files)} genomes")

    allg = []
    for f in files:
        iso = os.path.basename(f).replace(".tsv", "")
        g = per_gene(read_ips(f))
        g["isolate"] = iso
        allg.append(g)
        print(f"  {iso}: {len(g)} proteins with an InterPro entry, "
              f"{(g.go!='').sum()} with GO")
    A = pd.concat(allg, ignore_index=True)
    A.to_csv(f"{DATA}/interpro_annotation.csv", index=False)

    ref = A[A.isolate == REF].set_index("protein")
    if not len(ref):
        raise SystemExit(f"{REF} not annotated yet")

    # full proteome + lengths (genes with no InterPro hit must still count as background)
    L = {r.id: len(str(r.seq).rstrip("*"))
         for r in SeqIO.parse(f"analysis/mining/results/{REF}/{REF}.proteins.fa", "fasta")}
    p = pd.DataFrame({"protein": list(L), "plen": list(L.values())})
    p = p.merge(ref[["ipr", "ipr_desc", "go"]], on="protein", how="left").fillna("")
    print(f"\n{REF}: {len(p)} proteins; {(p.ipr!='').sum()} with InterPro "
          f"({100*(p.ipr!='').sum()/len(p):.1f}%), {(p.go!='').sum()} with GO "
          f"({100*(p.go!='').sum()/len(p):.1f}%)")
    print("  (PFAM-only annotation reached 67.3%)")

    v = pd.read_csv(f"{DATA}/clonal_variants{SFX}.csv").fillna("")
    imp = v[v.effect.isin(["missense", "nonsense", "stop_lost", "frameshift", "inframe_indel"])]
    hit = set(imp.transcript) - {""}
    p["hit"] = p.protein.isin(hit)
    print(f"variant-bearing genes: {p.hit.sum()}")

    rows = []
    for cat, terms in GO_CATS.items():
        p[cat] = p.go.map(lambda s: any(t in s.split("|") for t in terms))
    p["ANY_regulatory"] = p[list(GO_CATS)].any(axis=1)

    for cat in list(GO_CATS) + ["ANY_regulatory"]:
        m = p[cat]
        a, b = int((m & p.hit).sum()), int((~m & p.hit).sum())
        c, d = int((m & ~p.hit).sum()), int((~m & ~p.hit).sum())
        if a + c == 0:
            continue
        f = stats.fisher_exact([[a, b], [c, d]])
        # length-controlled logistic regression
        X = sm.add_constant(pd.DataFrame({"log_len": np.log10(p.plen),
                                          "cat": m.astype(int)}))
        try:
            mod = sm.Logit(p.hit.astype(int), X).fit(disp=0)
            p_adj, coef = mod.pvalues["cat"], mod.params["cat"]
        except Exception:
            p_adj, coef = np.nan, np.nan
        rows.append(dict(category=cat, n_hit=a, pct_hit=round(100*a/(a+b), 1),
                         n_bg=c, pct_bg=round(100*c/(c+d), 1),
                         odds=round(f.statistic, 2) if np.isfinite(f.statistic) else np.nan,
                         p_naive=round(f.pvalue, 4),
                         p_length_controlled=round(p_adj, 4), coef=round(coef, 3)))
    r = pd.DataFrame(rows)

    # length-matched permutation: random gene sets with the same protein-length decile
    # composition as the variant-bearing genes
    rng = np.random.default_rng(11)
    p["lbin"] = pd.qcut(p.plen, 10, labels=False, duplicates="drop")
    need = p[p.hit].lbin.value_counts()
    pools = {b: p.index[p.lbin == b].values for b in need.index}
    NPERM = 10000
    draws = [np.concatenate([rng.choice(pools[b], n, replace=False) for b, n in need.items()])
             for _ in range(NPERM)]
    r["p_length_perm"] = [
        round((sum(p.loc[d, c].sum() >= r_.n_hit for d in draws) + 1) / (NPERM + 1), 4)
        for c, r_ in zip(r.category, r.itertuples())]

    # biocontrol genes (CAZyme, MEROPS protease or SignalP-secreted) among variant-bearing genes
    a = pd.read_csv(f"analysis/mining/results/{REF}/{REF}.annotations.txt", sep="\t", dtype=str).fillna("")
    sec = set(open(f"analysis/mining/results/{REF}/{REF}.secreted_ids.txt").read().split())
    bc = set(a.TranscriptID[(a.CAZyme != "") | (a.Protease != "") | a.TranscriptID.isin(sec)])
    p["biocontrol"] = p.protein.isin(bc)
    k, n = int((p.biocontrol & p.hit).sum()), int(p.hit.sum())
    K, N = int(p.biocontrol.sum()), len(p)
    fb = stats.fisher_exact([[k, n - k], [K - k, N - n - (K - k)]])
    from statsmodels.stats.proportion import confint_proportions_2indep
    lo, hi = confint_proportions_2indep(k, n, K - k, N - n, method="newcomb")
    print(f"\nbiocontrol genes: {k}/{n} variant-bearing ({100*k/n:.1f}%) vs "
          f"{K-k}/{N-n} rest of proteome ({100*(K-k)/(N-n):.1f}%); "
          f"diff 95% CI {100*lo:+.1f} to {100*hi:+.1f} pp; Fisher P = {fb.pvalue:.3f}")
    r.to_csv(f"{DATA}/interpro_enrichment{SFX}.csv", index=False)
    print("\n=== GO-category enrichment among variant-bearing genes ===")
    print("(p_length_controlled = logistic regression with log10 protein length)")
    print(r.to_string(index=False))

    # annotated candidate list
    out = imp.merge(p[["protein", "plen", "ipr", "ipr_desc", "go", "ANY_regulatory"]],
                    left_on="transcript", right_on="protein", how="left")
    out.to_csv(f"{DATA}/clonal_variants_ipr{SFX}.csv", index=False)

    print("\n=== variant-bearing genes now characterised by InterPro "
          "(previously 'hypothetical protein') ===")
    named = out[(out.ipr_desc != "") & (out.ipr_desc.notna())].copy()
    named["ipr_desc"] = named.ipr_desc.str.slice(0, 78)
    print(named[["transcript", "effect", "aa_change", "ipr_desc"]]
          .drop_duplicates("transcript").to_string(index=False))
    print(f"\nWrote {DATA}/interpro_annotation.csv, clonal_variants_ipr.csv, "
          f"interpro_enrichment.csv")


if __name__ == "__main__":
    main()
