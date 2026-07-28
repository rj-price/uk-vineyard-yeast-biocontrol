import pandas as pd, glob, os
m=pd.read_csv("analysis/data/isolate_master.csv")
asco=m[m.species.isin(["Hanseniaspora uvarum","Pichia sp."])].copy()
fa={os.path.basename(f).split("_EKDN")[0]:f for f in glob.glob("genomics/*_scaffolds.fasta")}
gsize={"Hanseniaspora uvarum":9.0,"Pichia sp.":11.0}
rows=[]
for _,r in asco.iterrows():
    sid=r.seq_id
    old=os.path.abspath(fa[sid]) if sid in fa else "NA"
    rows.append((r.isolate,sid,r.species,gsize[r.species], old,
                 r.busco_complete, r.n50, r.size_inclN))
s=pd.DataFrame(rows,columns=["isolate","seq_id","species","exp_genome_mb",
                             "old_assembly","old_busco","old_n50","old_size"])
s=s.sort_values(["species","seq_id"]).reset_index(drop=True)
os.makedirs("analysis/reassembly",exist_ok=True)
s.to_csv("analysis/reassembly/samplesheet.tsv",sep="\t",index=False)
print(f"{len(s)} ascomycete isolates written")
print(s.to_string(index=False))
# pilot subset: one Hanseniaspora (Y2_16=RY3) + one Pichia (Y3_5=RY13)
pilot=s[s.seq_id.isin(["Y2_16","Y3_5"])]
pilot.to_csv("analysis/reassembly/samplesheet_pilot.tsv",sep="\t",index=False)
print("\nPilot:", list(pilot.seq_id))
