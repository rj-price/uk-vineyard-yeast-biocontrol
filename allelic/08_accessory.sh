#!/bin/bash
# Accessory gene presence/absence across all Kent H. uvarum (22) and Rhodotorula (8)
# (reviewer point 1). Only some isolates have funannotate gene models, so the
# annotated proteomes of each genus are pooled, clustered (CD-HIT 90%) into a pan
# set and projected onto every assembly with miniprot. Parsed by 36_accessory.py.
#SBATCH --job-name=accessory
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=0-1
#SBATCH --output=analysis/allelic/logs/accessory_%A_%a.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
E=/mnt/apps/users/jnprice/conda/envs
MP=$E/evm/bin/miniprot; CDHIT=$E/cdhit/bin/cd-hit
T=$SLURM_CPUS_PER_TASK
GRP=(Hanseniaspora Rhodotorula)
g=${GRP[$SLURM_ARRAY_TASK_ID]}
if [ "$g" = Hanseniaspora ]; then
  ANN=(RY3 RY4 RY7 RY14 RY25)
  ISO=(RY1 RY2 RY3 RY4 RY5 RY7 RY8 RY9 RY10 RY11 RY14 RY15 RY16 RY17 RY19 RY20 RY21 RY22 RY24 RY25 RY26 RY27)
else
  ANN=(UY1 TY27 TY28 TY29)
  ISO=(UY1 UY2 UY3 UY4 TY13 TY27 TY28 TY29)
fi
OUT=analysis/allelic/accessory/$g; mkdir -p "$OUT/paf"
: > "$OUT/pooled.faa"
for a in "${ANN[@]}"; do
  sed "s/^>/>${a}|/" analysis/mining/results/$a/$a.proteins.fa >> "$OUT/pooled.faa"
done
"$CDHIT" -i "$OUT/pooled.faa" -o "$OUT/pan.faa" -c 0.9 -aS 0.8 -n 5 -M 0 -T "$T" -d 0 > "$OUT/cdhit.log"
for i in "${ISO[@]}"; do
  "$MP" -t "$T" --outs=0.5 -N 1 ncbi_upload/genomes/$i.fasta "$OUT/pan.faa" > "$OUT/paf/$i.paf" 2> "$OUT/paf/$i.log"
  echo "$i $(cut -f1 "$OUT/paf/$i.paf" | sort -u | wc -l)"
done
echo DONE "$g"
