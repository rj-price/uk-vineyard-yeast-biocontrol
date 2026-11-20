#!/bin/bash
# All-vs-all fastANI on the eight Rhodotorula aff. babjevae isolates, which
# resolves the two near-clonal groups {UY1, UY4, TY27} and {TY28, UY2, UY3}.
#SBATCH --job-name=rhodo8ani
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=analysis/allelic/logs/rhodo8_ani.%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
L=analysis/allelic/Rhodotorula8_list.txt
for i in UY1 UY2 UY3 UY4 TY13 TY27 TY28 TY29; do
    echo "$PWD/ncbi_upload/genomes/$i.fasta"
done > "$L"
conda run -n fastani fastANI --ql "$L" --rl "$L" -t "$SLURM_CPUS_PER_TASK" \
    -o analysis/allelic/Rhodotorula8_ani.tsv
