#!/bin/bash
#SBATCH --job-name=ani2
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id/logs/ani2-%j.out
# Expanded whole-genome ANI: the 55 paper assemblies against the 30 reference
# genomes in species_id/refs/ (ref_final.tsv; Table S2), including the expanded
# Sporobolomyces panel. Output ani2_raw.tsv feeds Fig 2B (20_fig2_phylo.py),
# Fig S1 (28b_figS1_sporo_ani.py) and the TY7 vs V. victoriae value (85.4%).
# The fastANI call is the one recorded in SLURM accounting for job 30631357.
set -euo pipefail
SID=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id
cd "$SID"
mkdir -p logs
ls /mnt/shared/scratch/jnprice/private/biocontrol_yeast/ncbi_upload/genomes/*.fasta > ani2_queries.txt
ls refs/*.fna > ani2_refs.txt
conda run -n fastani fastANI --ql "$SID/ani2_queries.txt" --rl "$SID/ani2_refs.txt" \
    --matrix -t "${SLURM_CPUS_PER_TASK:-8}" -o "$SID/ani2_raw.tsv" 2> "$SID/logs/fastani2.log"
echo DONE
