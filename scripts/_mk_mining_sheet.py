import pandas as pd, glob, os
m=pd.read_csv("analysis/data/isolate_master.csv")
fa={os.path.basename(f).split("_EKDN")[0]:os.path.abspath(f) for f in glob.glob("genomics/*_scaffolds.fasta")}
plat={os.path.basename(f).split(".")[0]:os.path.abspath(f)
      for f in glob.glob("analysis/final_assemblies/ascomycetes_platanus/*.fasta")}

# one strong, phenotyped representative per species; assembly source chosen by quality
# (Metschnikowia assembled well with SPAdes; Hanseniaspora/Pichia use Platanus re-assemblies)
picks = [
    # isolate, seq_id, species_label, assembly_source
    ("UY1",  "Y1_1",  "Rhodotorula_babjevae",      "spades"),
    ("TY26", "Y2_9",  "Metschnikowia_ziziphicola", "spades"),
    ("TY30", "Y2_13", "Metschnikowia_pulcherrima", "spades"),
    ("TY4",  "Y1_8",  "Vishniacozyma_carnescens",  "spades"),
    ("TY6",  "Y1_10", "Sporobolomyces_roseus",     "spades"),
    ("RY14", "Y3_6",  "Hanseniaspora_uvarum",      "platanus"),
    ("RY13", "Y3_5",  "Pichia_sp",                 "platanus"),
]
rows=[]
for iso,sid,sp,src in picks:
    asm = fa[sid] if src=="spades" else plat[sid]
    assert os.path.exists(asm), f"missing {asm}"
    # AUGUSTUS pre-trained seed species (must be a valid `funannotate species` id),
    # chosen as the closest available relative for BUSCO-based augustus training:
    #   Metschnikowia -> candida_albicans (both Saccharomycotina CTG-clade)
    #   Hanseniaspora/Pichia -> saccharomyces_cerevisiae_S288C (Saccharomycotina)
    #   Vishniacozyma -> cryptococcus (both Tremellomycetes)
    #   Rhodotorula/Sporobolomyces -> ustilago_maydis (Basidiomycota)
    seed = {"Rhodotorula_babjevae":"ustilago_maydis",
            "Metschnikowia_ziziphicola":"candida_albicans",
            "Metschnikowia_pulcherrima":"candida_albicans",
            "Vishniacozyma_carnescens":"cryptococcus",
            "Sporobolomyces_roseus":"ustilago_maydis",
            "Hanseniaspora_uvarum":"saccharomyces_cerevisiae_S288C",
            "Pichia_sp":"saccharomyces_cerevisiae_S288C"}[sp]
    rows.append((iso,sid,sp,src,asm,seed))
s=pd.DataFrame(rows,columns=["isolate","seq_id","species","asm_source","assembly","seed_species"])
os.makedirs("analysis/mining",exist_ok=True)
s.to_csv("analysis/mining/samplesheet.tsv",sep="\t",index=False)
# pilot = the single most biocontrol-relevant genome (M. pulcherrima = pulcherrimin/iron)
s[s.isolate=="TY30"].to_csv("analysis/mining/samplesheet_pilot.tsv",sep="\t",index=False)
print(s.to_string(index=False))
print("\nPilot: TY30 (Metschnikowia pulcherrima)")
