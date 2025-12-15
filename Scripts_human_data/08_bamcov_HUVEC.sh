#!/bin/bash

set -euo pipefail

# --- Verbose flag (-v) ---
verbose=0
[[ "${1:-}" == "-v" ]] && verbose=1
vprint() { [[ $verbose -eq 1 ]] && echo "$@" || true; }

# --- Config ---
bam_dir="results/Human_cells/HUVEC_ATAC/mapping/"
bed_file="results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks/celltype_condition.filteredNfixed.union.peakSet"
UNION_PEAKSET="results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks/celltype_condition.filteredNfixed.union.peakSet"
UNION_SORTED="results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted.peakset"
output_file="results/Human_cells/HUVEC_ATAC/data/multicov_output_HUVEC.tsv"

# --- Checks & setup ---
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }
[[ -d "$bam_dir" ]] || { echo "ERROR: BAM dir not found: $bam_dir"; exit 1; }
[[ -f "$bed_file" ]] || { echo "ERROR: BED file not found: $bed_file"; exit 1; }
mkdir -p "$(dirname "$output_file")"

# Collect BAMs safely (handles 'no match' and spaces)
shopt -s nullglob
bam_files=( "${bam_dir}"*_rmdup.bam )
if (( ${#bam_files[@]} == 0 )); then
  echo "ERROR: No BAMs matching ${bam_dir}*_rmdup.bam"
  exit 1
fi


echo "🔧 Preparing union peakset..."
# Clean and sort union peakset
tail -n +2 "$UNION_PEAKSET" | awk '{OFS="\t"; print $1, $2, $3}' | sort -k1,1 -k2,2n > "$UNION_SORTED"


vprint "Starting bedtools multicov..."
vprint "BED file: $bed_file"
vprint "BAM files:"
for b in "${bam_files[@]}"; do vprint "  - $b"; done

# --- Run multicov ---
# Note: bedtools multicov prints BED columns + per-BAM counts in the order provided.
bedtools multicov -bams "${bam_files[@]}" -bed "$UNION_SORTED" > "$output_file"

vprint "Output written to $output_file"
vprint "Done."

