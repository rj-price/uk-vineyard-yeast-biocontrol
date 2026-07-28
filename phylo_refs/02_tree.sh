#!/bin/bash
# Build the reference-inclusive supermatrix from the BUSCO runs, then IQ-TREE
# with ultrafast bootstrap. Runs after the BUSCO array (afterok dependency).
#SBATCH --job-name=phylo_tree
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=analysis/phylo_refs/logs/tree_%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
DIR=analysis/phylo_refs
BP=/mnt/apps/users/jnprice/BUSCO_phylogenomics/BUSCO_phylogenomics.py

echo "=== completed BUSCO runs: $(ls -d $DIR/busco_runs/*/ | wc -l)"
rm -rf "$DIR/phylo_out"
# BUSCOs complete & single-copy in >=90% of taxa (100% too strict across 83 divergent taxa)
conda run -n BUSCO_phylogenomics python "$BP" \
    -i "$DIR/busco_runs" -o "$DIR/phylo_out" -t "$SLURM_CPUS_PER_TASK" \
    --supermatrix_only -psc 90

SM=$(find "$DIR/phylo_out" -name "SUPERMATRIX.fasta" | head -1)
echo "=== supermatrix: $SM"
grep -c '>' "$SM" || true

cd "$(dirname "$SM")"
conda run -n BUSCO_phylogenomics iqtree2 -s SUPERMATRIX.fasta \
    -m LG+F+G4 -B 1000 -T AUTO -ntmax "$SLURM_CPUS_PER_TASK" -pre refs_tree -redo
echo "DONE tree -> $(pwd)/refs_tree.contree"
