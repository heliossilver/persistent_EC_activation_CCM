#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
MERGED_PEAK_DIR="results/ATAC/merged_peaks"
OUTPUT_DIR="results/ATAC/filtered_peaksets"
UNION_PEAKSET="${MERGED_PEAK_DIR}/celltype_condition.filteredNfixed.union.peakSet"
UNION_SORTED="${MERGED_PEAK_DIR}/celltype_condition.filteredNfixed.union_clean_sorted.peakset"

FLOX_PEAK_DIR="${OUTPUT_DIR}/peakset_all_flox"
CELLTYPE_PEAK_DIR="${OUTPUT_DIR}/peakset_by_celltype"
FLOX_COMBINED_DIR="${OUTPUT_DIR}/combined_peak_all_flox"
CELLTYPE_COMBINED_DIR="${OUTPUT_DIR}/combined_peak_by_celltype"

CELLTYPES=("artery" "capillary_artery" "capillary" "capillary_vein" "large_vein")

mkdir -p "$FLOX_PEAK_DIR" "$CELLTYPE_PEAK_DIR" "$FLOX_COMBINED_DIR" "$CELLTYPE_COMBINED_DIR" "results/ATAC/data"

echo "🔧 Preparing union peakset..."
# Clean and sort union peakset
tail -n +2 "$UNION_PEAKSET" | awk '{OFS="\t"; print $1, $2, $3}' | sort -k1,1 -k2,2n > "$UNION_SORTED"

### STEP 1: All Flox Peakset
echo "📦 Creating combined Flox peakset..."
for celltype in "${CELLTYPES[@]}"; do
    tail -n +2 "${MERGED_PEAK_DIR}/${celltype}_Flox.filterNfixed.peakset"
done > "${FLOX_COMBINED_DIR}/all_flox_combined.peakset"

awk '{OFS="\t"; print $1, $2, $3}' "${FLOX_COMBINED_DIR}/all_flox_combined.peakset" | \
    sort -k1,1 -k2,2n > "${FLOX_COMBINED_DIR}/all_flox_combined_clean_sorted.peakset"

# Intersect with union
bedtools intersect -a "${FLOX_COMBINED_DIR}/all_flox_combined_clean_sorted.peakset" \
                   -b "$UNION_SORTED" -wb | \
    awk '{OFS="\t"; print $4, $5, $6}' | sort -k1,1 -k2,2n | uniq > "${FLOX_PEAK_DIR}/all_flox_filtered_union.peakset"

echo -e "chr\tstart\tend" | cat - "${FLOX_PEAK_DIR}/all_flox_filtered_union.peakset" > "${FLOX_PEAK_DIR}/all_flox_filtered_union.final.peakset"
rm "${FLOX_PEAK_DIR}/all_flox_filtered_union.peakset"

### STEP 2: By Cell Type (Flox + KO)
echo "🔁 Processing cell types individually..."
for celltype in "${CELLTYPES[@]}"; do
    echo "🧬 $celltype"

    tail -n +2 "${MERGED_PEAK_DIR}/${celltype}_Flox.filterNfixed.peakset" > tmp_flox
    tail -n +2 "${MERGED_PEAK_DIR}/${celltype}_KO.filterNfixed.peakset" > tmp_ko
    cat tmp_flox tmp_ko > "${CELLTYPE_COMBINED_DIR}/${celltype}_combined.peakset"

    awk '{OFS="\t"; print $1, $2, $3}' "${CELLTYPE_COMBINED_DIR}/${celltype}_combined.peakset" | \
        sort -k1,1 -k2,2n > "${CELLTYPE_COMBINED_DIR}/${celltype}_combined_clean_sorted.peakset"

    bedtools intersect -a "${CELLTYPE_COMBINED_DIR}/${celltype}_combined_clean_sorted.peakset" \
                       -b "$UNION_SORTED" -wb | \
        awk '{OFS="\t"; print $4, $5, $6}' | sort -k1,1 -k2,2n | uniq > "${CELLTYPE_PEAK_DIR}/${celltype}.filteredNfixed.union.peakset"

    echo -e "chr\tstart\tend" | cat - "${CELLTYPE_PEAK_DIR}/${celltype}.filteredNfixed.union.peakset" > \
        "${CELLTYPE_PEAK_DIR}/${celltype}.filteredNfixed.union.final.peakset"

    rm tmp_flox tmp_ko "${CELLTYPE_PEAK_DIR}/${celltype}.filteredNfixed.union.peakset"
done

echo "✅ All peaksets generated and filtered."


### STEP 3: Annotate union peakset
echo "🧬 Annotating union peakset..."

awk '{OFS="\t"; print $1, $2, $3, $1"-"$2"-"$3}' "$UNION_SORTED" > "${MERGED_PEAK_DIR}/celltype_condition.filteredNfixed.union_clean_sorted_label.peakset"

annotatePeaks.pl "${MERGED_PEAK_DIR}/celltype_condition.filteredNfixed.union_clean_sorted_label.peakset" mm10 -size 200 > results/ATAC/data/annotated_union_clean_sorted.peakset.txt

echo "📄 Annotation saved to: results/ATAC/data/annotated_union_clean_sorted.peakset.txt"



