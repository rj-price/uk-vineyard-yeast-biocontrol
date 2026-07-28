#!/bin/bash
#SBATCH --job-name=ani_within
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=analysis/allelic/logs/ani_%A_%a.log
#SBATCH --array=0-1
set -euo pipefail
GRP=(Rhodotorula Hanseniaspora)
g=${GRP[$SLURM_ARRAY_TASK_ID]}
conda run -n fastani fastANI --ql analysis/allelic/${g}_list.txt --rl analysis/allelic/${g}_list.txt \
  -t "$SLURM_CPUS_PER_TASK" -o analysis/allelic/${g}_ani_within.tsv
echo "DONE $g"
