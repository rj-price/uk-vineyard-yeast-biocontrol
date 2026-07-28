#!/bin/bash
#SBATCH --job-name=reasm_extract
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/reassembly/logs/%x-%j.out
set -euo pipefail
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
RA=$PROJ/analysis/reassembly
mkdir -p "$RA/reads"; cd "$RA/reads"
echo "[$(date)] extracting $(wc -l < "$RA/members.txt") read files (26 ascomycetes) ..."
tar -xvf "$PROJ/X204SC24101850-Z01-F001.tar" -T "$RA/members.txt" --strip-components=2
echo "[$(date)] done"; du -sh "$RA/reads"
