#!/bin/bash
#SBATCH --job-name=killer_tblastn
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/mining/logs/killer-%j.out
set -euo pipefail
# Targeted tblastn of a curated fungal killer-toxin panel against every assembly.
# PFAM missed these; killer toxins are often dsRNA/plasmid-borne or too divergent to
# sit in the predicted gene models, so we search the genomic DNA directly.
# One isolate-tagged nucleotide DB + a single tblastn (avoids 60x conda startup).
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
MINE=$PROJ/analysis/mining
QUERY=$MINE/killer_toxins.faa
SUBJ=$MINE/killer_subjects.tsv
OUT=$MINE/killer_toxin
T=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$OUT"

# 1. concatenate all assemblies, prefixing every contig with the isolate id
comb=$OUT/all_assemblies.fna
: > "$comb"
tail -n +2 "$SUBJ" | while IFS=$'\t' read -r iso sid sp asm; do
  [ -f "$asm" ] || { echo "MISSING $iso $asm" >&2; continue; }
  awk -v iso="$iso" '/^>/{sub(/^>/,">"iso"__");print;next}{print}' "$asm" >> "$comb"
done
echo "combined contigs: $(grep -c '>' "$comb")  from $(grep -c $'\t' "$SUBJ") subjects"

# 2. single tblastn against the isolate-tagged DB
conda run -n blast bash -c "
  makeblastdb -in '$comb' -dbtype nucl -out '$OUT/all_db' >/dev/null
  tblastn -query '$QUERY' -db '$OUT/all_db' -evalue 1e-4 -seg yes -num_threads $T \
    -outfmt '6 qseqid qlen sseqid pident length mismatch gapopen evalue bitscore qcovs' \
    > '$OUT/killer_tblastn.tsv'
"
echo "hit rows (e<1e-4): $(wc -l < "$OUT/killer_tblastn.tsv")"
echo "isolates with a hit: $(cut -f3 "$OUT/killer_tblastn.tsv" | sed 's/__.*//' | sort -u | wc -l)"
rm -f "$comb" "$OUT"/all_db.*   # tidy the big intermediate + db
