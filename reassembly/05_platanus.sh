#!/bin/bash
#SBATCH --job-name=platanus
#SBATCH --partition=medium
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=18:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/reassembly/logs/platanus-%A_%a.out
set -euo pipefail
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
RA=$PROJ/analysis/reassembly
SHEET=${1:-$RA/samplesheet.tsv}
TARGET_COV=${2:-70}
T=${SLURM_CPUS_PER_TASK:-24}; MEMGB=72

PLAT=$RA/bin/platanus_allee
FASTP=/mnt/apps/users/jnprice/conda/envs/polysnip-nf/bin/fastp
SEQTK=/mnt/apps/users/jnprice/conda/envs/ont_assembly/bin/seqtk
SEQKIT=/mnt/apps/users/jnprice/conda/envs/illumina_assembly/bin/seqkit
BUSCO_DL=/mnt/apps/users/jnprice/databases/busco_downloads

row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET")
isolate=$(echo "$row"|cut -f1); seq_id=$(echo "$row"|cut -f2)
species=$(echo "$row"|cut -f3); gmb=$(echo "$row"|cut -f4)
gbp=$(awk -v m="$gmb" 'BEGIN{printf "%d", m*1000000}')
echo "[$(date)] Platanus-allee $isolate ($seq_id) $species  ${TARGET_COV}x of ${gmb}Mb  T=$T"

OUT=$RA/platanus/$seq_id; mkdir -p "$OUT"
work=$RA/platanus/$seq_id/work; rm -rf "$work"; mkdir -p "$work"; cd "$work"
R1=$(ls "$RA/reads/$seq_id/"*_1.fq.gz); R2=$(ls "$RA/reads/$seq_id/"*_2.fq.gz)

# 1. trim
"$FASTP" -i "$R1" -I "$R2" -o t1.fq.gz -O t2.fq.gz \
  -j "$OUT/${seq_id}.fastp.json" -h "$OUT/${seq_id}.fastp.html" \
  -w "$T" --detect_adapter_for_pe 2> "$OUT/${seq_id}.fastp.log"

# 2. downsample to target coverage, output plain fastq (Platanus IP libs)
clean_bp=$("$SEQKIT" stats -T t1.fq.gz t2.fq.gz | awk 'NR>1{s+=$5} END{print s}')
frac=$(awk -v t="$TARGET_COV" -v g="$gbp" -v b="$clean_bp" 'BEGIN{f=t*g/b; print (f>=1)?1:f}')
echo "[$(date)] clean=$clean_bp (~$(awk -v g=$gbp -v b=$clean_bp 'BEGIN{printf "%.0f",b/g}')x) frac=$frac"
"$SEQTK" sample -s100 t1.fq.gz "$frac" > ds_1.fq
"$SEQTK" sample -s100 t2.fq.gz "$frac" > ds_2.fq
rm -f t1.fq.gz t2.fq.gz

# 3. Platanus-allee: assemble -> phase -> consensus
"$PLAT" assemble -f ds_1.fq ds_2.fq -t "$T" -m "$MEMGB" -o out -tmp "$work" 2> "$OUT/${seq_id}.assemble.log"
"$PLAT" phase -c out_contig.fa -IP1 ds_1.fq ds_2.fq -t "$T" -o out -tmp "$work" 2> "$OUT/${seq_id}.phase.log"
"$PLAT" consensus -c out_consensusInput.fa -IP1 ds_1.fq ds_2.fq -t "$T" -o out -tmp "$work" 2> "$OUT/${seq_id}.consensus.log"

# 4. final consensus haploid assembly
cons=$(ls out_consensusScaffold.fa 2>/dev/null || ls out_consensusScaffold*.fa | head -1)
"$SEQKIT" seq -m 500 "$cons" > "$OUT/${seq_id}.platanus.fasta"
"$SEQKIT" stats -a -T "$OUT/${seq_id}.platanus.fasta" > "$OUT/${seq_id}.stats.tsv"

# 5. BUSCO
conda run -n busco busco -i "$OUT/${seq_id}.platanus.fasta" -l fungi_odb10 \
  -o busco_${seq_id} -m genome -c "$T" -f --offline --download_path "$BUSCO_DL" \
  > "$OUT/${seq_id}.busco.log" 2>&1 || echo "BUSCO warn"
cp busco_${seq_id}/short_summary.*.txt "$OUT/" 2>/dev/null || true

cat "$OUT/${seq_id}.stats.tsv"; grep -hoE 'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%' "$OUT/"short_summary.*.txt 2>/dev/null | head -1
rm -rf "$work"      # drop reads + platanus intermediates, keep outputs
echo "[$(date)] done $seq_id"
