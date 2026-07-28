#!/bin/bash
#SBATCH --job-name=yeast_extract
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/mapping_qc/logs/%x-%j.out
set -euo pipefail

# Extract raw reads for the 15 target samples in a single pass through the tar.
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
QC=$PROJ/analysis/mapping_qc
TAR=$PROJ/X204SC24101850-Z01-F001.tar
READS=$QC/reads
mkdir -p "$READS"
cd "$READS"

echo "[$(date)] extracting $(wc -l < "$QC/members.txt") members from tar ..."
# GNU tar: -T lists exact member names to extract; strip the two leading dirs
tar -xvf "$TAR" -T "$QC/members.txt" \
    --strip-components=2
echo "[$(date)] extraction done"
ls -R "$READS" | head
du -sh "$READS"
