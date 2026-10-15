#!/bin/bash
# BUSCO (fungi_odb10, genome mode) on the 2 added Sporobolomyces references
# (pararoseus, salmonicolor) so they can enter the phylogenomic supermatrix tree
# alongside S. roseus. Outputs collect in busco_runs/ for BUSCO_phylogenomics.
#SBATCH --job-name=busco_sporo
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-2
#SBATCH --output=analysis/phylo_refs/logs/busco_sporo_%A_%a.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
DIR=analysis/phylo_refs
RUNS=$DIR/busco_runs
mkdir -p "$RUNS" "$DIR/logs"
DB=/mnt/apps/users/jnprice/databases/busco_downloads

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$DIR/genomes_sporo.tsv")   # +1 skips header
name=$(echo "$line" | cut -f1)
asm=$(echo "$line" | cut -f2)
[ -n "$name" ] || { echo "no task $SLURM_ARRAY_TASK_ID"; exit 0; }

if [ -f "$RUNS/$name/short_summary.specific.fungi_odb10.$name.txt" ]; then
    echo "$name already done"; exit 0
fi

echo "=== BUSCO $name : $asm"
conda run -n busco busco -i "$asm" -o "$name" --out_path "$RUNS" \
    -l fungi_odb10 -m genome -c "$SLURM_CPUS_PER_TASK" -f \
    --offline --download_path "$DB"
echo "DONE $name"
