#!/bin/bash
#SBATCH --job-name=reasm
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/reassembly/logs/reasm-%A_%a.out
set -euo pipefail

PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
RA=$PROJ/analysis/reassembly
SHEET=${1:-$RA/samplesheet.tsv}              # pass sheet as arg1 (avoids --export env issues)
TARGET_COV=${2:-70}                          # arg2
THREADS=${SLURM_CPUS_PER_TASK:-16}
MEMGB=48

FASTP=/mnt/apps/users/jnprice/conda/envs/polysnip-nf/bin/fastp
SEQTK=/mnt/apps/users/jnprice/conda/envs/ont_assembly/bin/seqtk
SEQKIT=/mnt/apps/users/jnprice/conda/envs/illumina_assembly/bin/seqkit
SPADES=/mnt/apps/users/jnprice/conda/envs/illumina_assembly/bin/spades.py
BUSCO_ENV=/mnt/apps/users/jnprice/conda/envs/busco
BUSCO_DL=/mnt/apps/users/jnprice/databases/busco_downloads

row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET")
isolate=$(echo "$row"|cut -f1); seq_id=$(echo "$row"|cut -f2)
species=$(echo "$row"|cut -f3); gmb=$(echo "$row"|cut -f4)
gbp=$(awk -v m="$gmb" 'BEGIN{printf "%d", m*1000000}')
echo "[$(date)] $isolate ($seq_id) $species  target ${TARGET_COV}x of ${gmb}Mb"

OUT=$RA/assemblies/$seq_id; mkdir -p "$OUT"
work=$(mktemp -d "${TMPDIR:-/tmp}/reasm_${seq_id}.XXXXXX"); trap 'rm -rf "$work"' EXIT
R1=$(ls "$RA/reads/$seq_id/"*_1.fq.gz); R2=$(ls "$RA/reads/$seq_id/"*_2.fq.gz)

# 1. adapter/quality trim
"$FASTP" -i "$R1" -I "$R2" -o "$work/t1.fq.gz" -O "$work/t2.fq.gz" \
   -j "$OUT/${seq_id}.fastp.json" -h "$OUT/${seq_id}.fastp.html" \
   -w "$THREADS" --detect_adapter_for_pe 2> "$OUT/${seq_id}.fastp.log"

# 2. downsample to TARGET_COV using post-trim base count
clean_bp=$("$SEQKIT" stats -T "$work/t1.fq.gz" "$work/t2.fq.gz" | awk 'NR>1{s+=$5} END{print s}')
frac=$(awk -v t="$TARGET_COV" -v g="$gbp" -v b="$clean_bp" 'BEGIN{f=t*g/b; print (f>=1)?1:f}')
raw_cov=$(awk -v g="$gbp" -v b="$clean_bp" 'BEGIN{printf "%.0f", b/g}')
echo "[$(date)] clean bases=$clean_bp (~${raw_cov}x) -> downsample fraction=$frac"
if awk -v f="$frac" 'BEGIN{exit !(f>=1)}'; then
  ds1="$work/t1.fq.gz"; ds2="$work/t2.fq.gz"; echo "already <= target, using all reads"
else
  "$SEQTK" sample -s100 "$work/t1.fq.gz" "$frac" | gzip > "$work/ds1.fq.gz"
  "$SEQTK" sample -s100 "$work/t2.fq.gz" "$frac" | gzip > "$work/ds2.fq.gz"
  ds1="$work/ds1.fq.gz"; ds2="$work/ds2.fq.gz"
fi

# 3. SPAdes (careful mode, multi-k) on node-local scratch
"$SPADES" --careful -k 21,33,55,77,99,127 -t "$THREADS" -m "$MEMGB" \
   -1 "$ds1" -2 "$ds2" -o "$work/spades" > "$OUT/${seq_id}.spades.log" 2>&1

# 4. keep contigs >=500 bp as the new assembly
"$SEQKIT" seq -m 500 "$work/spades/scaffolds.fasta" > "$OUT/${seq_id}.reassembled.fasta"
"$SEQKIT" stats -a -T "$OUT/${seq_id}.reassembled.fasta" > "$OUT/${seq_id}.stats.tsv"

# 5. BUSCO (offline, fungi_odb10) — must run inside the conda env
cd "$work"
conda run -n busco busco -i "$OUT/${seq_id}.reassembled.fasta" -l fungi_odb10 \
   -o busco_${seq_id} -m genome -c "$THREADS" -f --offline \
   --download_path "$BUSCO_DL" > "$OUT/${seq_id}.busco.log" 2>&1 || echo "BUSCO warn"
cp "$work/busco_${seq_id}/short_summary."*.txt "$OUT/" 2>/dev/null || true

echo "[$(date)] done $seq_id"
grep -E "num_seqs|sum_len|N50" "$OUT/${seq_id}.stats.tsv" || true
