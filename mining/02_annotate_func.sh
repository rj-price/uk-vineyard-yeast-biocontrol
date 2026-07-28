#!/bin/bash
#SBATCH --job-name=mine_annot
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis/mining/logs/annot-%A_%a.out
set -euo pipefail
PROJ=/mnt/shared/scratch/jnprice/private/biocontrol_yeast
MINE=$PROJ/analysis/mining
SHEET=${1:-$MINE/samplesheet.tsv}
T=${SLURM_CPUS_PER_TASK:-16}
export FUNANNOTATE_DB=/mnt/apps/users/jnprice/databases/funannotate_db
iso=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SHEET" | cut -f1)
predir=$MINE/work/$iso/predict
OUT=$MINE/results/$iso
cd "$MINE/work/$iso"
# funannotate annotate adds PFAM (incl. killer-toxin domains), MEROPS proteases,
# CAZymes & product names on top of the predict gene models.
# Skip GeneMark (expired) and eggNOG-mapper (needs ~44GB --dbmem; not needed for
# PFAM/MEROPS/CAZyme mining). shim/emapper.py fails the version check so funannotate
# treats eggNOG as not installed.
conda run -n funannotate bash -c "
  export GENEMARK_PATH=/tmp/no_genemark_annot
  export PATH=\$(echo \"\$PATH\" | tr ':' '\n' | grep -v gmes_linux | paste -sd:)
  funannotate annotate -i '$predir' --cpus $T --busco_db dikarya --eggnog $MINE/empty.emapper.annotations
" > "$OUT/07_annotate.log" 2>&1
# copy annotation table
cp "$predir"/../annotate_results/*.annotations.txt "$OUT/${iso}.annotations.txt" 2>/dev/null || \
  cp "$predir"/annotate_results/*.annotations.txt "$OUT/${iso}.annotations.txt" 2>/dev/null || \
  find "$MINE/work/$iso" -name "*.annotations.txt" -exec cp {} "$OUT/${iso}.annotations.txt" \; 2>/dev/null || true
echo "done $iso; annotations: $(wc -l < "$OUT/${iso}.annotations.txt" 2>/dev/null || echo NA)"
