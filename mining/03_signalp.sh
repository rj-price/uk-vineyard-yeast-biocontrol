#!/bin/bash
#SBATCH --job-name=signalp
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/mining/logs/signalp-%A_%a.out
set -euo pipefail
# Predict the SECRETED proteome (Sec/SPI signal peptides) for each mined genome
# with SignalP 5.0b, so the CAZyme/protease genotype can be restricted to the
# actually-secreted enzymes that the plate assays measure.
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
MINE=$PROJ/analysis/mining
SHEET=${1:-$MINE/samplesheet.tsv}
SIGNALP=/mnt/apps/users/jnprice/signalp/signalp-5.0b/bin/signalp
iso=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$SHEET" | cut -f1)  # +2: line 1 is header
RES=$MINE/results/$iso
faa=$RES/$iso.proteins.fa
out=$RES/signalp

mkdir -p "$out" "$out/tmp_$iso"
echo "[$iso] $(grep -c '>' "$faa") proteins -> SignalP 5.0b (euk, short)"
# short format = one summary table; -org euk for fungi
"$SIGNALP" -fasta "$faa" -org euk -format short -prefix "$out/$iso" -tmp "$out/tmp_$iso"
# extract transcript IDs predicted to carry a Sec/SPI signal peptide
awk -F'\t' '$0!~/^#/ && $2=="SP(Sec/SPI)"{print $1}' "$out/${iso}_summary.signalp5" > "$RES/${iso}.secreted_ids.txt"
echo "[$iso] secreted (SP) proteins: $(wc -l < "$RES/${iso}.secreted_ids.txt")"
