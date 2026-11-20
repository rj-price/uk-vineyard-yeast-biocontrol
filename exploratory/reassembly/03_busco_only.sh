#!/bin/bash
#SBATCH --job-name=reasm_busco
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/reassembly/logs/busco-%A_%a.out
set -euo pipefail
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
RA=$PROJ/analysis/reassembly
SHEET=${1:-$RA/samplesheet.tsv}
BUSCO_DL=/mnt/apps/users/jnprice/databases/busco_downloads
seq_id=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET" | cut -f2)
OUT=$RA/assemblies/$seq_id
work=$(mktemp -d "${TMPDIR:-/tmp}/busco_${seq_id}.XXXXXX"); trap 'rm -rf "$work"' EXIT
cd "$work"
conda run -n busco busco -i "$OUT/${seq_id}.reassembled.fasta" -l fungi_odb10 \
   -o busco_${seq_id} -m genome -c "${SLURM_CPUS_PER_TASK:-12}" -f --offline \
   --download_path "$BUSCO_DL" > "$OUT/${seq_id}.busco.log" 2>&1
cp "$work/busco_${seq_id}/short_summary."*.txt "$OUT/"
echo "done $seq_id"; grep -E "C:|complete" "$OUT/short_summary."*.txt | head -1
