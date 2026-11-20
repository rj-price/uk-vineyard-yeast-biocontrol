#!/bin/bash
# Re-run the clonal-pair variant call against the funannotate `sorted.fa` genomes
# so SNP coordinates are directly on the scaffolds the gene models use.
#SBATCH --job-name=clonal_ref
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=analysis/allelic/logs/clonalref_%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
OUT=analysis/allelic/clonal_ref; mkdir -p "$OUT"
for q in TY27 TY28 TY29; do
  echo "=== UY1 (ref, funannotate scaffolds) vs $q"
  conda run -n mummer dnadiff -p "$OUT/UY1_vs_$q" \
    analysis/mining/work/UY1/sorted.fa analysis/mining/work/$q/sorted.fa
  grep -E 'TotalSNPs|TotalIndels|AvgIdentity' "$OUT/UY1_vs_$q.report" | head -4
done
echo DONE
