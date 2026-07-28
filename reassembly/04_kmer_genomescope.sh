#!/bin/bash
#SBATCH --job-name=kmer_gs
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/reassembly/logs/kmer-%A_%a.out
set -euo pipefail
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
RA=$PROJ/analysis/reassembly
SHEET=${1:-$RA/samplesheet.tsv}
JF=/mnt/apps/users/jnprice/conda/envs/jellyfish/bin/jellyfish
GS=/mnt/apps/users/jnprice/conda/envs/jellyfish/bin/genomescope2
seq_id=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET" | cut -f2)
T=${SLURM_CPUS_PER_TASK:-12}
OUT=$RA/kmer/$seq_id; mkdir -p "$OUT"
work=$(mktemp -d "${TMPDIR:-/tmp}/kmer_${seq_id}.XXXXXX"); trap 'rm -rf "$work"' EXIT
R1=$(ls "$RA/reads/$seq_id/"*_1.fq.gz); R2=$(ls "$RA/reads/$seq_id/"*_2.fq.gz)

echo "[$(date)] jellyfish count k=21 $seq_id"
"$JF" count -C -m 21 -s 2G -t "$T" -o "$work/k21.jf" \
   <(zcat "$R1") <(zcat "$R2")
"$JF" histo -t "$T" "$work/k21.jf" > "$OUT/${seq_id}.histo"
echo "[$(date)] genomescope2"
"$GS" -i "$OUT/${seq_id}.histo" -o "$OUT" -k 21 -n "$seq_id" > "$OUT/${seq_id}.gs.log" 2>&1 || echo "GS warn"
echo "== ${seq_id} GenomeScope summary =="
cat "$OUT/${seq_id}_summary.txt" 2>/dev/null || cat "$OUT/summary.txt" 2>/dev/null || tail -20 "$OUT/${seq_id}.gs.log"
