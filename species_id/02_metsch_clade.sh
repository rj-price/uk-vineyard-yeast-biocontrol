#!/bin/bash
# fastANI of the two Metschnikowia isolates (TY26, TY30) against the public
# M. pulcherrima-clade genomes listed in metsch_clade/acc.txt (including the type
# material of M. fructicola and M. andauensis). Downloads with NCBI `datasets`.
#SBATCH --job-name=metschani
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=analysis/species_id/metsch_clade/ani.%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
D=analysis/species_id/metsch_clade
analysis/species_id/datasets download genome accession --inputfile "$D/acc.txt" \
    --include genome --filename "$D/dl.zip"
analysis/species_id/datasets summary genome accession --inputfile "$D/acc.txt" \
    --as-json-lines > "$D/meta.jsonl"
unzip -oq "$D/dl.zip" -d "$D/dl"
ls "$PWD/$D"/dl/ncbi_dataset/data/*/*.fna > "$D/refs.txt"
printf "%s\n" "$PWD/ncbi_upload/genomes/TY26.fasta" "$PWD/ncbi_upload/genomes/TY30.fasta" > "$D/q.txt"
conda run -n fastani fastANI --ql "$D/q.txt" --rl "$D/refs.txt" -t "$SLURM_CPUS_PER_TASK" \
    -o "$D/ani.tsv"
