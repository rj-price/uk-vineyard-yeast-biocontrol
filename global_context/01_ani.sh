#!/bin/bash
#SBATCH --job-name=gc_ani
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/global_context/logs/ani-%j.out
set -euo pipefail
GC=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/global_context
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$GC/logs"; cd "$GC"
for tag in rhodo huvarum; do
  echo "=== fastANI all-vs-all: $tag ==="
  conda run -n fastani fastANI --ql ${tag}_list.txt --rl ${tag}_list.txt \
    -t "$T" -o ${tag}_ani.tsv 2> logs/fastani_${tag}.log
  echo "$tag rows: $(wc -l < ${tag}_ani.tsv)"
done
echo "DONE"
