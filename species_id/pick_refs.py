#!/usr/bin/env python
"""Query NCBI for one best reference assembly per candidate species and emit an
accession list. Preference: type material > RefSeq > assembly level > N50."""
import json, subprocess, sys
DS="/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id/datasets"
# candidate species to place our isolates against (include both members of each
# discriminating pair, plus congeners that could be the true ID)
SPECIES=[
 "Hanseniaspora uvarum","Hanseniaspora opuntiae","Hanseniaspora guilliermondii",
 "Pichia kluyveri","Pichia fermentans","Pichia membranifaciens",
 "Metschnikowia pulcherrima","Metschnikowia ziziphicola","Metschnikowia fructicola","Metschnikowia chrysoperlae",
 "Rhodotorula babjevae","Rhodotorula graminis","Rhodotorula mucilaginosa","Rhodotorula toruloides",
 "Sporobolomyces roseus","Sporobolomyces japonicus",
 "Vishniacozyma carnescens","Vishniacozyma victoriae",
 "Filobasidium stepposum","Filobasidium magnum",
 "Hannaella oryzae","Hannaella sinensis","Hannaella zeae",
 "Naganishia vaughanmartiniae","Naganishia albida","Naganishia diffluens",
]
LEVEL={"Complete Genome":4,"Chromosome":3,"Scaffold":2,"Contig":1}

def summ(sp):
    try:
        out=subprocess.run([DS,"summary","genome","taxon",sp],capture_output=True,text=True,timeout=90).stdout
        return json.loads(out).get("reports",[]) or []
    except Exception as e:
        print("  ERR",sp,e,file=sys.stderr); return []

rows=[]
for sp in SPECIES:
    reps=summ(sp)
    if not reps:
        print(f"{sp}\tNONE",file=sys.stderr); continue
    def score(r):
        ai=r.get("assembly_info",{})
        typ=1 if ai.get("biosample",{}).get("attributes") else 0
        # type material flag
        froms=str(ai.get("assembly_info",{}))
        istype=1 if (ai.get("type_material") or "type" in str(r.get("type_material","")).lower()) else 0
        refseq=1 if r.get("accession","").startswith("GCF_") else 0
        lvl=LEVEL.get(ai.get("assembly_level",""),0)
        n50=r.get("assembly_stats",{}).get("contig_n50",0) or 0
        return (istype,refseq,lvl,n50)
    best=max(reps,key=score)
    ai=best.get("assembly_info",{})
    rows.append((sp,best["accession"],ai.get("assembly_level",""),
                 str(best.get("organism",{}).get("organism_name","")),
                 "TYPE" if ai.get("type_material") else ""))
with open("ref_accessions.tsv","w") as f:
    f.write("query_species\taccession\tlevel\torganism\ttype\n")
    for r in rows: f.write("\t".join(map(str,r))+"\n")
for r in rows: print("\t".join(map(str,r)))
print(f"\n{len(rows)} references selected",file=sys.stderr)
