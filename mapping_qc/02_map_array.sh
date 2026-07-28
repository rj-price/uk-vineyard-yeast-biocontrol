#!/bin/bash
#SBATCH --job-name=yeast_mapqc
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-15%5
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/mapping_qc/logs/mapqc-%A_%a.out
set -euo pipefail

PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
QC=$PROJ/analysis/mapping_qc
SHEET=$QC/samplesheet.tsv
READS=$QC/reads
RES=$QC/results
THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$RES" "$QC/logs"

FASTP=/mnt/apps/users/jnprice/conda/envs/polysnip-nf/bin/fastp
BWA=/mnt/apps/users/jnprice/conda/envs/variant_calling/bin/bwa
SAM=/mnt/apps/users/jnprice/conda/envs/variant_calling/bin/samtools

# --- pull this task's row (data line = array idx, +1 to skip header) ---
row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET")
isolate=$(echo "$row" | cut -f1); seq_id=$(echo "$row" | cut -f2)
species=$(echo "$row" | cut -f3); group=$(echo "$row" | cut -f4)
assembly=$(echo "$row" | cut -f6)
echo "[$(date)] $isolate ($seq_id) $species [$group]"

work=$(mktemp -d "${TMPDIR:-/tmp}/mapqc_${seq_id}.XXXXXX")
trap 'rm -rf "$work"' EXIT

R1=$(ls "$READS/$seq_id/"*_1.fq.gz)
R2=$(ls "$READS/$seq_id/"*_2.fq.gz)

# --- 1. fastp: adapter/quality trim + QC (duplication, adapter, Q30, GC) ---
"$FASTP" -i "$R1" -I "$R2" -o "$work/t1.fq.gz" -O "$work/t2.fq.gz" \
   -j "$RES/${seq_id}.fastp.json" -h "$RES/${seq_id}.fastp.html" \
   -w "$THREADS" --detect_adapter_for_pe 2> "$RES/${seq_id}.fastp.log"

# --- 2. map trimmed reads to the isolate's own assembly ---
cp "$assembly" "$work/asm.fa"
"$BWA" index "$work/asm.fa" 2>> "$RES/${seq_id}.bwa.log"
"$BWA" mem -t "$THREADS" "$work/asm.fa" "$work/t1.fq.gz" "$work/t2.fq.gz" 2>> "$RES/${seq_id}.bwa.log" \
  | "$SAM" sort -@ "$THREADS" -o "$work/aln.bam" -
"$SAM" index "$work/aln.bam"

# --- 3. mapping stats + per-contig coverage ---
"$SAM" flagstat "$work/aln.bam" > "$RES/${seq_id}.flagstat.txt"
"$SAM" coverage "$work/aln.bam" > "$RES/${seq_id}.coverage.txt"

echo "[$(date)] done $seq_id"
# work dir (reads-derived tmp + bam) auto-removed by trap
