#!/bin/bash
# BUSCO (fungi_odb10, genome mode) on 83 genomes = 60 isolates + 21 references +
# 2 Fusarium outgroups, so the phylogenomic tree can include the ANI reference
# genomes shown in Fig 2B. Outputs collect in busco_runs/ for BUSCO_phylogenomics.
#SBATCH --job-name=busco_refs
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=1-83%6
#SBATCH --output=analysis/phylo_refs/logs/busco_%A_%a.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
DIR=analysis/phylo_refs
RUNS=$DIR/busco_runs
mkdir -p "$RUNS"
DB=/mnt/apps/users/jnprice/databases/busco_downloads

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$DIR/genomes.tsv")   # +1 skips header
name=$(echo "$line" | cut -f1)
asm=$(echo "$line" | cut -f2)
[ -n "$name" ] || { echo "no task $SLURM_ARRAY_TASK_ID"; exit 0; }

# resume-safe: skip if this genome's BUSCO already completed
if [ -f "$RUNS/$name/short_summary.specific.fungi_odb10.$name.txt" ]; then
    echo "$name already done"; exit 0
fi

echo "=== BUSCO $name : $asm"
conda run -n busco busco -i "$asm" -o "$name" --out_path "$RUNS" \
    -l fungi_odb10 -m genome -c "$SLURM_CPUS_PER_TASK" -f \
    --offline --download_path "$DB"
echo "DONE $name"
