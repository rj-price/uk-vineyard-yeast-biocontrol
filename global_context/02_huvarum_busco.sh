#!/bin/bash
# BUSCO (fungi_odb10, genome mode) on the 19 global H. uvarum reference genomes.
# Our 22 Kent H. uvarum isolates already have BUSCO runs (phylo_refs/busco_runs);
# they are symlinked in at the tree step, so only the refs need computing here.
#SBATCH --job-name=huv_busco
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-19%6
#SBATCH --output=analysis/global_context/logs/huvbusco_%A_%a.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
GC=analysis/global_context
RUNS=$GC/huvarum_busco
mkdir -p "$RUNS" "$GC/logs"
DB=/mnt/apps/users/jnprice/databases/busco_downloads

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$GC/huvarum_refs_manifest.tsv")
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
