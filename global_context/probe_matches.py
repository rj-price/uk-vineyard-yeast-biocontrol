import json, pandas as pd

d = json.load(open("analysis/global_context/meta/rhodotorula_genus.json"))
want = {"GCA_026119225.1", "GCA_920103735.2", "GCA_937859235.1",
        "GCA_920104245.2", "GCA_001329695.1"}  # 3 conspecific + babjevae T + graminis WP1


def attrs(rep):
    bs = rep.get("assembly_info", {}).get("biosample", {})
    return {x.get("name"): x.get("value") for x in bs.get("attributes", [])}


for r in d["reports"]:
    if r["accession"] not in want:
        continue
    org = r.get("organism", {}); ai = r.get("assembly_info", {}); a = attrs(r)
    bs = ai.get("biosample", {})
    print("=" * 70)
    print(r["accession"], "|", org.get("organism_name"),
          "| strain:", (org.get("infraspecific_names", {}) or {}).get("strain"))
    print("  submitter :", ai.get("submitter"))
    print("  biosample :", bs.get("description", {}).get("title"))
    for k in ["geo_loc_name", "country", "isolation_source", "isolation source",
              "host", "collection_date", "collected_by", "type-material"]:
        if a.get(k):
            print(f"  {k}: {a[k]}")
    print("  bioproject:", [p.get("accession") for p in
                            (bs.get("bioprojects") or [])])
