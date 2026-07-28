#!/bin/bash
# Assemble a combined BUSCO dir (22 Kent isolate runs symlinked from phylo_refs +
# 19 global ref runs) and build the H. uvarum supermatrix tree with IQ-TREE UFBoot.
#SBATCH --job-name=huv_tree
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=analysis/global_context/logs/huvtree_%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
GC=analysis/global_context
BP=/mnt/apps/users/jnprice/BUSCO_phylogenomics/BUSCO_phylogenomics.py
COMB=$GC/huvarum_busco_all
rm -rf "$COMB"; mkdir -p "$COMB"

# 19 global refs (computed here)
for d in "$GC"/huvarum_busco/*/; do
    [ -d "$d" ] && ln -s "$(readlink -f "$d")" "$COMB/$(basename "$d")"
done
# 22 Kent isolates (reuse existing runs)
while read -r name; do
    [ -n "$name" ] || continue
    src="analysis/phylo_refs/busco_runs/$name"
    [ -d "$src" ] && ln -s "$(readlink -f "$src")" "$COMB/$name"
done < "$GC/huvarum_isolates.txt"

echo "=== combined BUSCO taxa: $(ls -d $COMB/*/ | wc -l)"
rm -rf "$GC/huvarum_phylo"
conda run -n BUSCO_phylogenomics python "$BP" \
    -i "$COMB" -o "$GC/huvarum_phylo" -t "$SLURM_CPUS_PER_TASK" \
    --supermatrix_only -psc 90

SM=$(find "$GC/huvarum_phylo" -name "SUPERMATRIX.fasta" | head -1)
echo "=== supermatrix: $SM  ($(grep -c '>' "$SM") taxa)"
cd "$(dirname "$SM")"
conda run -n BUSCO_phylogenomics iqtree2 -s SUPERMATRIX.fasta \
    -m LG+F+G4 -B 1000 -T AUTO -ntmax "$SLURM_CPUS_PER_TASK" -pre huvarum_tree -redo
echo "DONE tree -> $(pwd)/huvarum_tree.contree"
