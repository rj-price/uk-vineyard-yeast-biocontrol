#!/bin/bash
# Direct variant call between the clonal pair UY1 (Bc+) and TY27 (Bc-), 99.99% ANI.
# At this divergence the genome-wide variant list is small enough to enumerate
# exhaustively, so we do not need an outlier test - every difference is a candidate.
# TY28/TY29 (~98.75% from UY1) are included as outgroups: a variant that is
# UY1-specific vs ALL THREE is more interesting than one shared with TY28/29.
#SBATCH --job-name=clonal_pair
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=analysis/allelic/logs/clonal_%j.log
set -euo pipefail

PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
cd "$PROJ"
OUT=analysis/allelic/clonal
mkdir -p "$OUT"

get_asm () { awk -F'\t' -v i="$1" '$1==i{print $4}' analysis/mining/killer_subjects.tsv; }
UY1=$(get_asm UY1)

for q in TY27 TY28 TY29; do
  asm=$(get_asm "$q")
  echo "=== UY1 (ref) vs $q"
  conda run -n mummer dnadiff -p "$OUT/UY1_vs_$q" "$UY1" "$asm"
  echo "--- $q report:"
  grep -E 'TotalSNPs|TotalIndels|AlignedBases|AvgIdentity' "$OUT/UY1_vs_$q.report" || true
done

echo DONE
