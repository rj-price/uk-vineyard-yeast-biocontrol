#!/bin/bash
#SBATCH --job-name=species_ani
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id/logs/spid-%j.out
set -euo pipefail
SID=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$SID/logs"; cd "$SID"
# whole-genome ANI: 60 isolate assemblies vs 19 species reference genomes
conda run -n fastani fastANI --ql query_list.txt --rl ref_list.txt \
  --matrix -t "$T" -o ani_raw.tsv 2> logs/fastani.log
echo "ANI rows: $(wc -l < ani_raw.tsv)"
echo "DONE"
