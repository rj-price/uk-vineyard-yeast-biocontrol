#!/bin/bash
# Callable genome for the clonal variant calls (07_clonal_call.sh + 38_clonal_reads.py):
# all-sites haploid calls, counting positions where every member is called with
# DP >= 10 and >= 90% of reads on the called allele (same filters as 38, minus QUAL,
# which is not defined for invariant sites).
#SBATCH --job-name=callable
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=analysis/allelic/logs/callable_%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
OUT=analysis/allelic/clonal_reads
BCF=/mnt/apps/users/jnprice/conda/envs/variant_calling/bin/bcftools
count() {  # ref, members...
  local ref=$1; shift
  local bams=(); for s in "$@"; do bams+=("$OUT/${s}_on_${ref}.bam"); done
  local glen; glen=$(awk '!/^>/{n+=length($0)} END{print n}' "$OUT/ref/$ref.fa")
  "$BCF" mpileup --threads 4 -f "$OUT/ref/$ref.fa" -q 20 -Q 20 -a FORMAT/AD -Ou "${bams[@]}" \
    | "$BCF" call --threads 4 --ploidy 1 -m -Ou \
    | "$BCF" query -f '[%GT:%AD\t]\n' \
    | awk -v ref="$ref" -v glen="$glen" -F'\t' '
        { ok=1
          for (i=1; i<NF; i++) {
            split($i, a, ":"); if (a[1]=="." ) {ok=0; break}
            n=split(a[2], ad, ","); dp=0; for (j=1;j<=n;j++) dp+=ad[j]
            if (dp < 10 || ad[a[1]+1] < 0.9*dp) {ok=0; break}
          }
          c+=ok }
        END { printf "%s\t%d\t%d\t%.1f\n", ref, c, glen, 100*c/glen }'
}
printf "group\tcallable_bp\tassembly_bp\tpct\n" > analysis/data/clonal_callable.tsv
count UY1 UY1 TY27 UY4 >> analysis/data/clonal_callable.tsv
count TY28 TY28 UY2 UY3 >> analysis/data/clonal_callable.tsv
cat analysis/data/clonal_callable.tsv
