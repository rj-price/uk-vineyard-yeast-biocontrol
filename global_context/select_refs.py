import json, pandas as pd

META = "analysis/global_context/meta"


def st(r): return r.get("assembly_stats", {})


def attrs(rep):
    bs = rep.get("assembly_info", {}).get("biosample", {})
    return {x.get("name"): x.get("value") for x in bs.get("attributes", [])}


def geo_of(a):
    return (a.get("geo_loc_name") or a.get("country")
            or a.get("geographic location (country and/or sea)"))


def rowify(r):
    org = r.get("organism", {}); ai = r.get("assembly_info", {}); a = attrs(r)
    return dict(acc=r["accession"], sp=org.get("organism_name"),
                strain=(org.get("infraspecific_names", {}) or {}).get("strain"),
                level=ai.get("assembly_level"), n50=st(r).get("contig_n50") or 0,
                geo=geo_of(a), iso_src=a.get("isolation_source") or a.get("isolation source"),
                host=a.get("host"), date=a.get("collection_date") or a.get("collection date"))


# ---- Rhodotorula curated set ----
d = json.load(open(f"{META}/rhodotorula_genus.json"))
df = pd.DataFrame([rowify(r) for r in d["reports"]])
df["species"] = df.sp.str.extract(r"^Rhodotorula (\w+)")[0]
# "sp" = unnamed R. sp. strains — keep ALL (a cryptic lineage could match ours)
take_all = {"babjevae", "graminis", "glutinis", "diobovata", "sphaerocarpa",
            "araucariae", "dairenensis", "taiwanensis", "sp"}
sel = []
for sp, grp in df.groupby("species", dropna=True):
    grp = grp.sort_values("n50", ascending=False)
    if sp in take_all:
        sel.append(grp)
    else:                            # distant relatives / large sets -> 2 best each
        sel.append(grp.head(2))
rh = pd.concat(sel).drop_duplicates("acc")
rh.to_csv(f"{META}/rhodotorula_refs.csv", index=False)
print("Rhodotorula refs:", len(rh))
print(rh.sp.value_counts().to_string())

# ---- H. uvarum full set ----
dh = json.load(open(f"{META}/huvarum.json"))
hf = pd.DataFrame([rowify(r) for r in dh["reports"]])
hf.to_csv(f"{META}/huvarum_refs.csv", index=False)
print("\nH. uvarum refs:", len(hf))
