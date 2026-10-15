#!/bin/bash
# Joint haploid variant calling per clonal group on the BAMs from 06_clonal_reads.sh.
#SBATCH --job-name=clonal_call
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=analysis/allelic/logs/clonalcall_%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
OUT=analysis/allelic/clonal_reads
BCF=/mnt/apps/users/jnprice/conda/envs/variant_calling/bin/bcftools
call() {  # ref, members...
  local ref=$1; shift
  local bams=(); for s in "$@"; do bams+=("$OUT/${s}_on_${ref}.bam"); done
  "$BCF" mpileup --threads 8 -f "$OUT/ref/$ref.fa" -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -Ou "${bams[@]}" \
    | "$BCF" call --threads 8 --ploidy 1 -mv -Oz -o "$OUT/group_$ref.vcf.gz"
  "$BCF" index -f "$OUT/group_$ref.vcf.gz"
}
call UY1 UY1 TY27 UY4
call TY28 TY28 UY2 UY3
echo DONE
