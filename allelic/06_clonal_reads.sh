#!/bin/bash
# Read-level validation of within-clone variants (reviewer point 3).
# Two near-clonal Rhodotorula groups (fastANI ~99.99%):
#   A: UY1 (ref), TY27, UY4      -> reads mapped to UY1 funannotate sorted.fa
#   B: TY28 (ref), UY2, UY3      -> reads mapped to TY28 funannotate sorted.fa
# Every isolate, including the reference itself, is mapped, so self-calls in the
# reference sample give the assembly-error baseline.
# Array index 0-5 = one isolate; after the array, run 07_clonal_call.sh.
#SBATCH --job-name=clonal_reads
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-5
#SBATCH --output=analysis/allelic/logs/clonalreads_%A_%a.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
ISO=(UY1 TY27 UY4 TY28 UY2 UY3)
SEQ=(Y1_1 Y2_10 Y1_4 Y2_11 Y1_2 Y1_3)
REF=(UY1 UY1 UY1 TY28 TY28 TY28)
i=$SLURM_ARRAY_TASK_ID; iso=${ISO[$i]}; sid=${SEQ[$i]}; ref=${REF[$i]}
OUT=analysis/allelic/clonal_reads; mkdir -p "$OUT/ref"
E=/mnt/apps/users/jnprice/conda/envs
FASTP=$E/polysnip-nf/bin/fastp; BWA=$E/variant_calling/bin/bwa; SAM=$E/variant_calling/bin/samtools
T=$SLURM_CPUS_PER_TASK

work=$(mktemp -d "${TMPDIR:-/tmp}/clr_${iso}.XXXXXX"); trap 'rm -rf "$work"' EXIT
tar -xf X204SC24101850-Z01-F001.tar -C "$work" --wildcards "X204SC24101850-Z01-F001/01.RawData/${sid}/${sid}_*.fq.gz"
R1=$(ls "$work"/X204SC24101850-Z01-F001/01.RawData/$sid/*_1.fq.gz)
R2=$(ls "$work"/X204SC24101850-Z01-F001/01.RawData/$sid/*_2.fq.gz)

# reference index (each ref built once, guarded by a lock)
fa=$OUT/ref/$ref.fa
( flock 9
  if [ ! -s "$fa.bwt" ]; then
    cp analysis/mining/work/$ref/sorted.fa "$fa"
    "$BWA" index "$fa"; "$SAM" faidx "$fa"
  fi ) 9>"$OUT/ref/$ref.lock"

"$FASTP" -i "$R1" -I "$R2" -o "$work/t1.fq.gz" -O "$work/t2.fq.gz" -w "$T" \
  --detect_adapter_for_pe -j "$OUT/$iso.fastp.json" -h /dev/null 2> "$OUT/$iso.fastp.log"
"$BWA" mem -t "$T" -R "@RG\tID:$iso\tSM:$iso\tPL:ILLUMINA" "$fa" "$work/t1.fq.gz" "$work/t2.fq.gz" \
  | "$SAM" fixmate -m -@ "$T" - - \
  | "$SAM" sort -@ "$T" -o "$work/s.bam" -
"$SAM" markdup -@ "$T" "$work/s.bam" "$OUT/${iso}_on_${ref}.bam"
"$SAM" index "$OUT/${iso}_on_${ref}.bam"
"$SAM" flagstat "$OUT/${iso}_on_${ref}.bam" > "$OUT/${iso}_on_${ref}.flagstat"
echo DONE "$iso"
