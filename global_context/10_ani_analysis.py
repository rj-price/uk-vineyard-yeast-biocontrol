"""Parse the fastANI all-vs-all output for each genus and answer:
  Rhodotorula — what is the closest *described* species to the Kent lineage, and does
                anything (incl. unnamed R. sp.) reach the 95% species threshold?
  H. uvarum   — confirm the Kent isolates are conspecific with global H. uvarum, and
                summarise Kent-vs-Kent / Kent-vs-global ANI.
Writes tidy tables under analysis/global_context/data/.
"""
import os, pandas as pd, numpy as np

GC = "analysis/global_context"
os.makedirs(f"{GC}/data", exist_ok=True)


def load(tag):
    names = pd.read_csv(f"{GC}/{tag}_names.csv")
    p2n = dict(zip(names.path.map(os.path.abspath), names.name))
    p2k = dict(zip(names.path.map(os.path.abspath), names.kind))
    a = pd.read_csv(f"{GC}/{tag}_ani.tsv", sep="\t", header=None,
                    names=["q", "r", "ani", "frag", "tot"])
    a["qn"] = a.q.map(lambda p: p2n.get(os.path.abspath(p)))
    a["rn"] = a.r.map(lambda p: p2n.get(os.path.abspath(p)))
    a["qk"] = a.q.map(lambda p: p2k.get(os.path.abspath(p)))
    a["rk"] = a.r.map(lambda p: p2k.get(os.path.abspath(p)))
    return a, names


def species_of(refname, annot):
    row = annot[annot.acc == refname]
    if len(row):
        return str(row.sp.iloc[0])
    return refname


# ---------- Rhodotorula ----------
def rhodo():
    a, names = load("rhodo")
    annot = pd.read_csv(f"{GC}/meta/rhodotorula_annot.csv")
    _summary(a, "rhodo")
    # best reference hit per Kent isolate
    ref = a[(a.qk == "isolate") & (a.rk == "reference")].copy()
    ref["ref_sp"] = ref.rn.map(lambda n: species_of(n, annot))
    best = (ref.sort_values("ani", ascending=False)
            .groupby("qn").head(1)[["qn", "rn", "ref_sp", "ani"]]
            .rename(columns={"qn": "isolate", "rn": "best_ref", "ani": "best_ani"}))
    best = best.sort_values("best_ani", ascending=False)
    best.to_csv(f"{GC}/data/rhodo_best_hit.csv", index=False)
    print("\n===== Rhodotorula: closest described relative per Kent isolate =====")
    print(best.to_string(index=False))
    print(f"\nmax ANI of any Kent isolate to any reference: {ref.ani.max():.2f}%")
    print(f"  -> {'>=95 (conspecific!)' if ref.ani.max()>=95 else '<95 => NO described species matches: distinct/undescribed lineage'}")
    # Kent internal cohesion
    kk = a[(a.qk == "isolate") & (a.rk == "isolate") & (a.qn != a.rn)]
    print(f"Kent-vs-Kent ANI: {kk.ani.min():.2f}-{kk.ani.max():.2f}% "
          f"(median {kk.ani.median():.2f})  -> single tight lineage" )
    # per-species best (does ANY species, incl unnamed R. sp., reach 95?)
    ref["described"] = ~ref.ref_sp.str.contains(r"\bsp\.", regex=True)
    sp_best = ref.groupby("ref_sp").ani.max().sort_values(ascending=False)
    print("\nbest Kent->reference ANI per reference species (top 12):")
    print(sp_best.head(12).round(2).to_string())
    sp_best.round(3).to_csv(f"{GC}/data/rhodo_species_bestANI.csv")


# ---------- H. uvarum ----------
def huvarum():
    a, names = load("huvarum")
    kk = a[(a.qk == "isolate") & (a.rk == "isolate") & (a.qn != a.rn)]
    kr = a[(a.qk == "isolate") & (a.rk == "reference")]
    print("\n===== H. uvarum =====")
    print(f"Kent-vs-Kent ANI:   {kk.ani.min():.2f}-{kk.ani.max():.2f}% (median {kk.ani.median():.2f})")
    print(f"Kent-vs-global ANI: {kr.ani.min():.2f}-{kr.ani.max():.2f}% (median {kr.ani.median():.2f})")
    print(f"  -> all {'>=95 conspecific with global H. uvarum' if kr.ani.min()>=95 else 'CHECK: some <95'}")
    kr.groupby("qn").ani.max().round(3).to_csv(f"{GC}/data/huvarum_kent_bestglobal.csv")


def _summary(a, tag):
    n_pairs = len(a)
    print(f"[{tag}] fastANI reported {n_pairs} pairs; "
          f"{a.qn.isna().sum()} unmapped query, {a.rn.isna().sum()} unmapped ref")


if __name__ == "__main__":
    rhodo()
    huvarum()
