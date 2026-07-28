#!/bin/bash
#SBATCH --job-name=of_allelic
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=analysis/allelic/logs/of_%A_%a.log
#SBATCH --array=0-1
set -euo pipefail
GRP=(Rhodotorula Hanseniaspora)
g=${GRP[$SLURM_ARRAY_TASK_ID]}
d=analysis/allelic/$g
rm -rf "$d/OrthoFinder"
conda run -n orthofinder orthofinder -f "$d" -t "$SLURM_CPUS_PER_TASK" -a 2 -og
echo "DONE $g"
