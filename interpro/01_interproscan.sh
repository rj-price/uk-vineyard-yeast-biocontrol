#!/bin/bash
# Full functional annotation (InterPro + GO) of all 12 mined genomes.
# The funannotate annotations carry PFAM only - the InterPro/GO/EggNog columns are
# empty because the original functional-annotation step never delivered.
#
# NOTE: everything must live under the project dir. The session scratchpad is on
# node-local /tmp and is NOT visible to compute nodes.
#SBATCH --job-name=iprscan
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
# benchmark: 580s / 200 proteins on 16 cpu -> ~5h for a ~6k-protein genome (the
# benchmark's JVM+DB startup is a fixed cost, so the real figure should be lower).
# %3 caps this at 48 cpu / 96G concurrently.
#SBATCH --time=12:00:00
#SBATCH --array=0-11%3
#SBATCH --output=analysis/interpro/logs/ips_%A_%a.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
IPS=/mnt/apps/users/jnprice/interproscan-5.74-105.0/interproscan.sh
IPRO=analysis/interpro

iso=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$IPRO/isolates.txt")
[ -n "$iso" ] || { echo "no isolate for task $SLURM_ARRAY_TASK_ID"; exit 1; }
out="$IPRO/results/${iso}.tsv"

if [ -s "$out" ]; then echo "$iso already done, skipping"; exit 0; fi

tmp="$IPRO/tmp/$iso"; rm -rf "$tmp"; mkdir -p "$tmp"
echo "=== $iso : $(grep -c '>' $IPRO/${iso}.faa) proteins, $SLURM_CPUS_PER_TASK cpu"
S=$(date +%s)
# -dp: compute locally rather than querying the EBI lookup service (no external dependency)
conda run -n interproscan "$IPS" \
  -i "$IPRO/${iso}.faa" -f TSV -o "$out" \
  -goterms -iprlookup -dp -cpu "$SLURM_CPUS_PER_TASK" -t p \
  -T "$tmp"
echo "DONE $iso in $(( $(date +%s) - S ))s -> $(wc -l < $out) rows"
rm -rf "$tmp"
