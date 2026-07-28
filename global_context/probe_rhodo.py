import pandas as pd, os
GC = "analysis/global_context"
names = pd.read_csv(f"{GC}/rhodo_names.csv")
p2n = dict(zip(names.path.map(os.path.abspath), names.name))
a = pd.read_csv(f"{GC}/rhodo_ani.tsv", sep="\t", header=None,
                names=["q", "r", "ani", "frag", "tot"])
a["qn"] = a.q.map(lambda p: p2n.get(os.path.abspath(p)))
a["rn"] = a.r.map(lambda p: p2n.get(os.path.abspath(p)))
annot = pd.read_csv(f"{GC}/meta/rhodotorula_annot.csv")
sp = dict(zip(annot.acc, annot.sp)); geo = dict(zip(annot.acc, annot.country))
strain = dict(zip(annot.acc, annot.strain)); n50 = dict(zip(annot.acc, annot.n50))

ref = a[(a.qn == "UY1") & (a.rn.astype(str).str.startswith("GCA"))].copy()
ref["sp"] = ref.rn.map(sp); ref["country"] = ref.rn.map(geo)
ref["strain"] = ref.rn.map(strain); ref["n50"] = ref.rn.map(n50)
print("=== UY1 vs all Rhodotorula refs (ANI desc, top 20) ===")
print(ref.sort_values("ani", ascending=False)
      [["rn", "sp", "strain", "country", "n50", "ani", "frag", "tot"]]
      .head(20).to_string(index=False))

print("\n=== babjevae + graminis genomes: best ANI from ANY Kent isolate ===")
kent = a[a.qn.isin(names[names.kind == "isolate"].name)]
for acc in annot[annot.sp.str.contains("babjevae|graminis", case=False)].acc:
    sub = kent[kent.rn == acc]
    b = f"{sub.ani.max():.2f}" if len(sub) else "no hit (<~80%)"
    print(f"{acc}  {str(sp[acc]):26} strain={str(strain.get(acc)):12} "
          f"{str(geo.get(acc)):10} best={b}")
