"""Harmonise country + isolation-source metadata for the global reference genomes,
so tree tips / ANI rows can be annotated by geography and habitat."""
import pandas as pd, re

META = "analysis/global_context/meta"


def country(geo):
    if not isinstance(geo, str) or not geo.strip() or geo.lower() in (
            "none", "missing", "not applicable", "not collected"):
        return "unknown"
    return geo.split(":")[0].strip()


def source_cat(row):
    txt = " ".join(str(row.get(k, "")) for k in ("iso_src", "host", "strain")).lower()
    if any(w in txt for w in ["grape", "vitis", "wine", "must", "vineyard", "muscatel",
                              "chardonnay", "fermentation", "cherry"]):
        return "grape / wine / fruit ferment"
    if any(w in txt for w in ["kimchi", "food", "citrus", "kombucha"]):
        return "other food / plant"
    if any(w in txt for w in ["drosophila", "cepaea", "snail", "excrement", "insect",
                              "gut", "frass"]):
        return "insect / animal"
    if any(w in txt for w in ["soil", "leaf", "phyllo", "water", "marine", "sea",
                              "bark", "environment", "air"]):
        return "environmental"
    return "unknown / other"


for tag in ["rhodotorula", "huvarum"]:
    df = pd.read_csv(f"{META}/{tag}_refs.csv")
    df["country"] = df.geo.map(country)
    df["source"] = df.apply(source_cat, axis=1)
    df.to_csv(f"{META}/{tag}_annot.csv", index=False)
    print(f"===== {tag} ({len(df)}) =====")
    print("countries:", df.country.value_counts().to_dict())
    print("sources  :", df.source.value_counts().to_dict())
    print()
