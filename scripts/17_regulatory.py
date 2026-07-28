#!/usr/bin/env python
"""
17_regulatory.py
The clonal pair UY1 (Bc+) / TY27 (Bc-) differ by ~600 SNPs, with no enrichment of
variants in biocontrol genes. Test the alternative hypothesis: are the
protein-altering variants in TRANS-ACTING REGULATORS (transcription factors,
chromatin modifiers, signalling), which would change expression of a shared,
conserved biocontrol gene set?

Uses PFAM domains (67% of genes annotated) mapped to names via the funannotate
Pfam-A.clans.tsv. NOTE: InterPro/GO/EggNog columns are EMPTY in these
annotations (the functional-annotation step did not deliver), so PFAM is the
only domain evidence available.

-> analysis/data/clonal_regulatory.csv
"""
import os, re
import pandas as pd, numpy as np
from scipy import stats

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast"
DATA = f"{PROJ}/analysis/data"
PFAM_CLANS = "/mnt/apps/users/jnprice/databases/funannotate_db/Pfam-A.clans.tsv"
REF = "UY1"

# PFAM name/description patterns for trans-acting regulation, by category.
# Deliberately broad: we are testing a category-level hypothesis, not one gene.
CATS = {
    "TF_DNA_binding": [
        r"\bZn_clus\b", r"fungal.*transcription", r"Fungal_trans", r"\bGATA\b",
        r"\bbZIP", r"Homeobox", r"Myb_DNA-binding", r"\bzf-C2H2", r"HSF_DNA-bind",
        r"helix-turn-helix", r"\bHTH\b", r"Transcription factor", r"TF_",
        r"DNA-binding", r"\bAPSES\b", r"\bCP2\b", r"\bE2F", r"\bSTE\b",
        r"Copper-fist", r"\bZn\(2\)-C6\b",
    ],
    "chromatin_epigenetic": [
        r"Histone", r"\bSET\b", r"Chromo", r"Bromodomain", r"\bJmjC\b", r"\bJmjN\b",
        r"SNF2", r"Sir2", r"\bHDAC\b", r"Acetyltransferase.*histone", r"MOZ_SAS",
        r"\bPHD\b", r"Tudor", r"\bDNMT\b", r"DNA_methylase", r"Methyltransf.*histone",
        r"\bRRM\b.*chromatin", r"CAF1", r"ISWI", r"\bHMG\b",
    ],
    "signalling": [
        r"Pkinase", r"Protein kinase", r"\bPP2C\b", r"Phosphatase", r"Response_reg",
        r"\bGPCR\b", r"\bRGS\b", r"\bRas\b", r"Ras family", r"\bRho\b", r"\bGTPase\b",
        r"\bcNMP\b", r"Adenylate.*cyclase", r"\bPAS\b", r"His_Kinase", r"HATPase",
        r"\b14-3-3\b",
    ],
    "RNA_post_transcriptional": [
        r"\bRRM_", r"RNA recognition", r"\bKH_\b", r"PUF\b", r"Pumilio",
        r"\bLSM\b", r"DEAD", r"Helicase", r"CCCH", r"Poly\(A\)",
    ],
}


def load_pfam_names():
    """Pfam-A.clans.tsv: acc, clan, clan_name, short_name, description."""
    d = pd.read_csv(PFAM_CLANS, sep="\t", header=None, dtype=str,
                    names=["acc", "clan", "clan_name", "short", "desc"]).fillna("")
    d["label"] = d.short + " | " + d.desc
    return dict(zip(d.acc, d.label))


def classify(label):
    hits = []
    for cat, pats in CATS.items():
        if any(re.search(p, label, re.I) for p in pats):
            hits.append(cat)
    return hits


def main():
    os.chdir(PROJ)
    names = load_pfam_names()
    ann = pd.read_csv(f"analysis/mining/results/{REF}/{REF}.annotations.txt",
                      sep="\t", dtype=str).fillna("")

    def domains(pf):
        accs = [a.split(".")[0] for a in re.split(r"[;,\s]+", pf) if a.startswith("PF")]
        return [names.get(a, a) for a in accs]

    ann["dom"] = ann.PFAM.map(domains)
    ann["dom_str"] = ann.dom.map(lambda L: "; ".join(L))
    ann["cats"] = ann.dom_str.map(lambda s: classify(s) if s else [])
    ann["is_reg"] = ann.cats.map(len) > 0

    v = pd.read_csv(f"{DATA}/clonal_variants.csv").fillna("")
    imp = v[v.effect.isin(["missense", "nonsense", "stop_lost", "indel_in_gene"])]
    hit_tx = set(imp.transcript) - {""}

    ann["hit"] = ann.TranscriptID.isin(hit_tx)
    print(f"genes with protein-altering variants (UY1 vs TY27): {ann.hit.sum()}")
    print(f"PFAM-annotated: {(ann.PFAM!='').sum()}/{len(ann)} "
          f"({100*(ann.PFAM!='').sum()/len(ann):.1f}%)\n")

    # enrichment: regulatory genes among variant-bearing vs the rest of the proteome
    print("=== enrichment of regulatory categories among variant-bearing genes ===")
    rows = []
    for cat in list(CATS) + ["ANY_regulatory"]:
        if cat == "ANY_regulatory":
            mask = ann.is_reg
        else:
            mask = ann.cats.map(lambda L: cat in L)
        a = int((mask & ann.hit).sum())        # regulatory & hit
        b = int((~mask & ann.hit).sum())       # non-reg & hit
        c = int((mask & ~ann.hit).sum())       # regulatory & not hit
        d = int((~mask & ~ann.hit).sum())
        f = stats.fisher_exact([[a, b], [c, d]])
        rows.append(dict(category=cat, n_hit_genes=a, n_hit_total=a + b,
                         pct_hit=round(100 * a / (a + b), 1),
                         n_bg=c, pct_bg=round(100 * c / (c + d), 1),
                         odds=round(f.statistic, 2) if np.isfinite(f.statistic) else np.nan,
                         p=round(f.pvalue, 4)))
    r = pd.DataFrame(rows)
    print(r.to_string(index=False))

    # the actual regulatory genes carrying variants
    cand = ann[ann.hit & ann.is_reg].copy()
    out = imp.merge(cand[["TranscriptID", "dom_str", "cats"]],
                    left_on="transcript", right_on="TranscriptID", how="inner")
    out.to_csv(f"{DATA}/clonal_regulatory.csv", index=False)
    print(f"\n=== variant-bearing genes with regulatory domains: {len(cand)} ===")
    if len(cand):
        show = out[["transcript", "effect", "aa_change", "cats", "dom_str"]].copy()
        show["dom_str"] = show.dom_str.str.slice(0, 70)
        print(show.to_string(index=False))

    print(f"\n=== ALL variant-bearing genes with ANY PFAM domain "
          f"({(ann.hit & (ann.PFAM!='')).sum()}/{ann.hit.sum()}) ===")
    hp = ann[ann.hit & (ann.PFAM != "")][["TranscriptID", "dom_str"]]
    hp = hp.merge(imp[["transcript", "effect", "aa_change"]],
                  left_on="TranscriptID", right_on="transcript")
    hp["dom_str"] = hp.dom_str.str.slice(0, 62)
    print(hp[["TranscriptID", "effect", "aa_change", "dom_str"]].to_string(index=False))
    print(f"\nWrote {DATA}/clonal_regulatory.csv")


if __name__ == "__main__":
    main()
