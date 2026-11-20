#!/bin/bash
#SBATCH --job-name=d1d2
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/species_id/logs/d1d2-%j.out
# LSU D1/D2 barcode identification for isolates with no conspecific reference genome
# (Vishniacozyma TY4, TY7, TY10; Filobasidium TY14; Rhodotorula UY1 and TY13 as
# checks). Results feed scripts/14_species_final.py and Table S5 ("D1/D2 rDNA").
#
#   1. barrnap (eukaryotic models) extracts rRNA genes from each assembly; --outseq
#      writes them strand-corrected.
#   2. The longest 28S copy per isolate is kept (rdna/lsu_all.fna, <isolate>_28S).
#   3. Its first 650 bp, the 5' end of the LSU that spans the D1/D2 domains, is the
#      barcode (rdna/d1d2.fna).
#   4. blastn against NCBI nt (remote), top 5 hits per barcode (rdna/d1d2_blast.tsv).
#      The full-length LSU was also searched (rdna/lsu_blast.tsv) as a check.
#
# Run interactively in July 2026 and written up here afterwards. The commands
# reproduce the sequences in rdna/lsu_all.fna and rdna/d1d2.fna exactly (checked
# Nov 2026, from the NCBI-submitted copies of the assemblies); remote BLAST hits depend on
# the state of nt at the time of the search.
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
SID=analysis/species_id
OUT=$SID/rdna
mkdir -p "$OUT" "$SID/logs"
ISOLATES=(UY1 TY13 TY4 TY7 TY10 TY14)
T=${SLURM_CPUS_PER_TASK:-4}

: > "$OUT/lsu_all.fna"
for i in "${ISOLATES[@]}"; do
    conda run -n extract_its barrnap --kingdom euk --threads "$T" --quiet \
        --outseq "$OUT/${i}_rrna.fna" "ncbi_upload/genomes/$i.fasta" > /dev/null
    # longest 28S record
    awk -v id="${i}_28S" '
        /^>/ {keep = ($0 ~ /^>28S_rRNA/); if (keep) {n++}; next}
        keep {s[n] = s[n] $0}
        END {b = 1; for (k in s) if (length(s[k]) > length(s[b])) b = k
             print ">" id; print s[b]}' "$OUT/${i}_rrna.fna" >> "$OUT/lsu_all.fna"
done

# D1/D2 = first 650 bp of the 28S
awk '/^>/{print; next} {print substr($0, 1, 650)}' "$OUT/lsu_all.fna" > "$OUT/d1d2.fna"

FMT="6 qseqid pident length evalue bitscore stitle"
conda run -n blast blastn -remote -db nt -query "$OUT/d1d2.fna" \
    -max_target_seqs 5 -outfmt "$FMT" > "$OUT/d1d2_blast.tsv" 2> "$OUT/d1d2_blast.err"
conda run -n blast blastn -remote -db nt -query "$OUT/lsu_all.fna" \
    -outfmt "$FMT" > "$OUT/lsu_blast.tsv" 2> "$OUT/lsu_blast.err"
echo DONE
