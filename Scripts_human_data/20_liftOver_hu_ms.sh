#!/bin/bash
set -euo pipefail

# Directories and filenames
HUMAN_PEAKS="results/Human_D3/bulk_ATAC/merged_peaks/condition_treatment.filteredNfixed.union_clean_sorted_label.peakset"
MOUSE_PEAKS="results/ATAC/merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted_label.peakset"
CHAIN_HG38_MM10="reference/hg38ToMm10.over.chain"
CHAIN_MM10_HG38="reference/mm10ToHg38.over.chain"


# Output directories
OUTDIR="results/human_mouse_integration/liftover"
mkdir -p "$OUTDIR"

echo "▶️ Liftover human peaks to mm10..."
liftOver "$HUMAN_PEAKS" \
         "$CHAIN_HG38_MM10" \
         "$OUTDIR/human_peaks_mm10coord.bed" \
         "$OUTDIR/unmapped_hu_ms.bed" \
         -minMatch=0.5

echo "▶️ Liftover mouse peaks to hg38..."
liftOver "$MOUSE_PEAKS" \
         "$CHAIN_MM10_HG38" \
         "$OUTDIR/mouse_peaks_hg38coord.bed" \
         "$OUTDIR/unmapped_ms_hu.bed" \
         -minMatch=0.5

echo "🔎 Filtering lifted peaks for length ≤ 1000 and standard chromosomes..."
awk '($3 - $2 <= 1000) && ($1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/)' \
    "$OUTDIR/human_peaks_mm10coord.bed" > "$OUTDIR/human_peaks_mm10coord.filtered.bed"

awk '($3 - $2 <= 1000) && ($1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/)' \
    "$OUTDIR/mouse_peaks_hg38coord.bed" > "$OUTDIR/mouse_peaks_hg38coord.filtered.bed"

echo "🔁 Reciprocal liftover back to original genome..."
liftOver "$OUTDIR/human_peaks_mm10coord.filtered.bed" \
         "$CHAIN_MM10_HG38" \
         "$OUTDIR/human_peaks_backtohg38.bed" \
         "$OUTDIR/unmapped_back_hu.bed" \
         -minMatch=0.5

liftOver "$OUTDIR/mouse_peaks_hg38coord.filtered.bed" \
         "$CHAIN_HG38_MM10" \
         "$OUTDIR/mouse_peaks_backtomm10.bed" \
         "$OUTDIR/unmapped_back_ms.bed" \
         -minMatch=0.5

echo "🧹 Filtering for reciprocal peak ID match and standard chromosomes..."
awk 'NR==FNR {ids[$4]; next} $4 in ids && $1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/' \
    "$OUTDIR/human_peaks_backtohg38.bed" "$OUTDIR/human_peaks_mm10coord.filtered.bed" > "$OUTDIR/human_peaks_mm10coord.reciprocal.bed"

awk 'NR==FNR {ids[$4]; next} $4 in ids && $1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/' \
    "$OUTDIR/mouse_peaks_backtomm10.bed" "$OUTDIR/mouse_peaks_hg38coord.filtered.bed" > "$OUTDIR/mouse_peaks_hg38coord.reciprocal.bed"

echo "🔗 Intersecting with actual peak calls in other species..."
bedtools intersect -wo \
  -a "$OUTDIR/human_peaks_mm10coord.reciprocal.bed" \
  -b "$MOUSE_PEAKS" > "$OUTDIR/human_peaks_with_mouse_ortho_peak.tsv"

bedtools intersect -wo \
  -a "$OUTDIR/mouse_peaks_hg38coord.reciprocal.bed" \
  -b "$HUMAN_PEAKS" > "$OUTDIR/mouse_peaks_with_human_ortho_peak.tsv"

echo "✅ Done! Orthologous peak matches saved in: $OUTDIR"

