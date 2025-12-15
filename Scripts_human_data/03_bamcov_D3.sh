#!/usr/bin/env bash
set -euo pipefail

# --- Verbose flag (-v) ---
verbose=0
[[ "${1:-}" == "-v" ]] && verbose=1
vprint() { [[ $verbose -eq 1 ]] && echo "$@" || true; }

# ------- CONFIG -------
OLD="results/Human_cells/D3_ATAC/mergedpeaks/condition_treatment.filteredNfixed.union_clean_sorted_label.peakset"
NEW="results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted_label.peakset"

MASTER_V2="results/Human_cells/D3_ATAC/D3_merged_peaks/master_v2_union.peakset"               
MASTER_V2_SORTED="results/Human_cells/D3_ATAC/D3_merged_peaks/master_v2_union_clean_sorted.peakset"
MASTER_V2_SORTED_LABEL="results/Human_cells/D3_ATAC/D3_merged_peaks/master_v2_union_clean_sorted_label.peakset"
ORDER_FILE="results/Human_cells/D3_ATAC/data/multicov_bam_order_master2.txt"

DO_COUNTS=1
BAM_DIR="results/Human_cells/D3_ATAC/mapping/"
OUT_COUNTS="results/Human_cells/D3_ATAC/data/multicov_master_v2.tsv"

# --- Checks & setup ---
[[ -f "$MASTER_V2" ]] || { echo "ERROR: MASTER_V2 not found: $MASTER_V2"; exit 1; }
mkdir -p "$(dirname "$MASTER_V2_SORTED")" "$(dirname "$MASTER_V2_SORTED_LABEL")" \
         "$(dirname "$ORDER_FILE")" "$(dirname "$OUT_COUNTS")"

# Clean and sort union peakset (header-safe)
tmp=$(mktemp)
awk 'BEGIN{OFS="\t"}
     NR==1 && $2 !~ /^[0-9]+$/ {next}     # drop header if start col not integer
     {print $1,$2,$3}' "$MASTER_V2" \
| sort -k1,1 -k2,2n > "$tmp"
mv "$tmp" "$MASTER_V2_SORTED"

# Add label column chr-start-end
awk '{OFS="\t"; print $1, $2, $3, $1"-"$2"-"$3}' "$MASTER_V2_SORTED" > "$MASTER_V2_SORTED_LABEL"

echo "✅ master_v2 sorted:       $MASTER_V2_SORTED"
echo "✅ master_v2 sorted+label: $MASTER_V2_SORTED_LABEL"

# ---  run multicov ---
if [[ "$DO_COUNTS" -eq 1 ]]; then
  command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }
  [[ -d "$BAM_DIR" ]] || { echo "ERROR: BAM dir not found: $BAM_DIR"; exit 1; }

  shopt -s nullglob
  bam_files=( "${BAM_DIR}"*_rmdup.bam )
  (( ${#bam_files[@]} > 0 )) || { echo "ERROR: No BAMs matching ${BAM_DIR}*_rmdup.bam"; exit 1; }

  vprint "Starting bedtools multicov..."
  vprint "BAM files:"; for b in "${bam_files[@]}"; do vprint "  - $b"; done

  printf "%s\n" "${bam_files[@]}" > "$ORDER_FILE"
  vprint "BAM order written to $ORDER_FILE"

  bedtools multicov -bams "${bam_files[@]}" -bed "$MASTER_V2_SORTED" > "$OUT_COUNTS"
  echo "✅ Counts written: $OUT_COUNTS"
fi

vprint "Done."

