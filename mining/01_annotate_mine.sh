#!/bin/bash
#SBATCH --job-name=mine
#SBATCH --partition=medium
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/mining/logs/mine-%A_%a.out
set -euo pipefail

PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
MINE=$PROJ/analysis/mining
SHEET=${1:-$MINE/samplesheet.tsv}
T=${SLURM_CPUS_PER_TASK:-24}
export FUNANNOTATE_DB=/mnt/apps/users/jnprice/databases/funannotate_db
# GeneMark deliberately not used (Apr-2025 key likely expired); funannotate falls
# back to AUGUSTUS+SNAP+GlimmerHMM+EVM, which is ample for CAZyme/cluster mining.
DBCAN_DB=/mnt/apps/users/jnprice/databases/dbcan_db

row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET")
iso=$(echo "$row"|cut -f1); sid=$(echo "$row"|cut -f2); sp=$(echo "$row"|cut -f3)
asm=$(echo "$row"|cut -f5); seed=$(echo "$row"|cut -f6)
tag=$(echo "$iso" | tr -d '_')                     # locus_tag prefix
OUT=$MINE/results/$iso; mkdir -p "$OUT" "$MINE/logs"
work=$MINE/work/$iso; rm -rf "$work"; mkdir -p "$work"; cd "$work"
echo "[$(date)] MINE $iso ($sid) $sp  seed=$seed  T=$T"

# ---------------------------------------------------------------------------
# 1. funannotate: clean -> sort -> mask  (short contig headers, softmask repeats)
# ---------------------------------------------------------------------------
conda run -n funannotate funannotate clean -i "$asm" -o clean.fa 2> "$OUT/01_clean.log" || cp "$asm" clean.fa
conda run -n funannotate funannotate sort -i clean.fa -o sorted.fa -b scaffold --minlen 200 2> "$OUT/02_sort.log"
conda run -n funannotate funannotate mask -i sorted.fa -o masked.fa --cpus "$T" 2> "$OUT/03_mask.log"

# ---------------------------------------------------------------------------
# 2. funannotate predict -> gene models + proteome + GBK
# ---------------------------------------------------------------------------
sp_name=$(echo "$sp" | tr '_' ' ')
# funannotate finds GeneMark via ~/.bashrc but its licence is EXPIRED; hide it so
# predict falls back to AUGUSTUS+SNAP+GlimmerHMM+EVM (ample for CAZyme/cluster mining).
conda run -n funannotate bash -c "
  export GENEMARK_PATH=/tmp/no_genemark_$sid
  export PATH=\$(echo \"\$PATH\" | tr ':' '\n' | grep -v gmes_linux | paste -sd:)
  funannotate predict -i masked.fa -o predict -s '$sp_name' \
    --busco_seed_species $seed --busco_db dikarya --name $tag \
    --cpus $T --tmpdir '$work/tmp'
" > "$OUT/04_predict.log" 2>&1

prot=$(ls predict/predict_results/*.proteins.fa | head -1)
gbk=$(ls predict/predict_results/*.gbk | head -1)
cp "$prot" "$OUT/${iso}.proteins.fa"
cp "$gbk"  "$OUT/${iso}.gbk"
cp predict/predict_results/*.gff3 "$OUT/${iso}.gff3" 2>/dev/null || true
nprot=$(grep -c "^>" "$OUT/${iso}.proteins.fa" || echo 0)
echo "[$(date)] predicted $nprot proteins"

# ---------------------------------------------------------------------------
# 2b. funannotate annotate -> PFAM (killer-toxin domains), MEROPS proteases,
#     CAZymes, product names. GeneMark + eggNOG(--dbmem ~44GB) hidden via shim.
# ---------------------------------------------------------------------------
conda run -n funannotate bash -c "
  export GENEMARK_PATH=/tmp/no_genemark_$sid
  export PATH=\$(echo \"\$PATH\" | tr ':' '\n' | grep -v gmes_linux | paste -sd:)
  funannotate annotate -i predict --cpus $T --busco_db dikarya --eggnog $MINE/empty.emapper.annotations
" > "$OUT/07_annotate.log" 2>&1 || echo "annotate warn (see log)"
adir=predict/annotate_results
cp $adir/*.annotations.txt "$OUT/${iso}.annotations.txt" 2>/dev/null || \
  find predict -name "*.annotations.txt" -exec cp {} "$OUT/${iso}.annotations.txt" \; 2>/dev/null || true
# prefer the functionally-annotated GBK/proteins for downstream steps
 agbk=$(ls $adir/*.gbk 2>/dev/null | head -1); [ -n "$agbk" ] && cp "$agbk" "$OUT/${iso}.gbk"
aprot=$(ls $adir/*.proteins.fa 2>/dev/null | head -1); [ -n "$aprot" ] && cp "$aprot" "$OUT/${iso}.proteins.fa"

# ---------------------------------------------------------------------------
# 3. dbCAN -> CAZymes (lytic enzymes: chitinase/glucanase/protease etc.)
# ---------------------------------------------------------------------------
conda run -n run_dbcan run_dbcan CAZyme_annotation \
  --input_raw_data "$OUT/${iso}.proteins.fa" --mode protein \
  --output_dir "$OUT/dbcan" --db_dir "$DBCAN_DB" --threads "$T" \
  > "$OUT/05_dbcan.log" 2>&1 || echo "dbcan warn (see log)"

# ---------------------------------------------------------------------------
# 4. antiSMASH -> secondary metabolite clusters (siderophores, pulcherrimin, terpenes)
# ---------------------------------------------------------------------------
conda run -n antismash antismash --taxon fungi --genefinding-tool none \
  --output-dir "$OUT/antismash" --cpus "$T" \
  --cb-general --cb-knownclusters --pfam2go "$OUT/${iso}.gbk" \
  > "$OUT/06_antismash.log" 2>&1 || echo "antismash warn (see log)"

echo "[$(date)] done $iso : $nprot proteins, dbCAN=$(ls -d "$OUT/dbcan" 2>/dev/null && echo ok), antiSMASH=$(ls "$OUT/antismash/index.html" 2>/dev/null && echo ok)"
rm -rf "$work/tmp"      # keep predict/ trainers? drop tmp only
