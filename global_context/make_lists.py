import pandas as pd, os

GC = "analysis/global_context"
m = pd.read_csv("analysis/data/isolate_master.csv")
g = pd.read_csv("analysis/phylo_refs/genomes.tsv", sep="\t")
path_of = dict(zip(g.name, g.assembly))

sets = {
    "rhodo": (m.species_genome.str.contains("babjevae", case=False, na=False),
              f"{GC}/refs_rhodo"),
    "huvarum": (m.species_genome.str.contains("Hanseniaspora", case=False, na=False),
                f"{GC}/refs_huvarum"),
}

for tag, (mask, refdir) in sets.items():
    rows = []
    listf = open(f"{GC}/{tag}_list.txt", "w")
    for iso in m[mask].isolate:
        p = path_of.get(iso)
        if p and os.path.exists(p):
            listf.write(p + "\n")
            rows.append((os.path.abspath(p), iso, "Kent, UK", "isolate"))
    for fn in sorted(os.listdir(refdir)):
        if fn.endswith(".fna"):
            p = os.path.abspath(f"{refdir}/{fn}")
            listf.write(p + "\n")
            rows.append((p, fn[:-4], "", "reference"))
    listf.close()
    pd.DataFrame(rows, columns=["path", "name", "geo", "kind"]).to_csv(
        f"{GC}/{tag}_names.csv", index=False)
    n_iso = sum(1 for r in rows if r[3] == "isolate")
    print(f"{tag}: {len(rows)} genomes ({n_iso} Kent isolates + {len(rows)-n_iso} refs)")
