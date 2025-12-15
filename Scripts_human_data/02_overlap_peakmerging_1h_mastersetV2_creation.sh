#!/bin/bash
set -euo pipefail  # safer shell

# === CONFIGURATION ===
BLACKLIST="reference/hg38-blacklist.v2.bed.gz"
CHROMSIZE="reference/hg38.chrom.sizes"
OUT_BASE="results/Human_cells/D3_ATAC"
MERGED_PEAK_DIR="${OUT_BASE}/D3_merged_peaks"
PEAK_LIST_FILE="${OUT_BASE}/folder_summit_paths.txt"
MACS2_OUT_BASE="${OUT_BASE}/macs2"
MERGING_SCRIPT="scripts/utils/iterative_overlap_peak_merging.R"
GENOME="hg38"
UNION_PEAKSET="results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union.peakSet"
UNION_SORTED="results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted.peakset"

# Create output folders
mkdir -p "$MERGED_PEAK_DIR"
mkdir -p "$(dirname "$PEAK_LIST_FILE")"
mkdir -p "results/logs"


# Generate list of .narrowPeak files
echo "🧾 Creating peak list file: $PEAK_LIST_FILE"
echo -n "" > "$PEAK_LIST_FILE"


for folder in "$MACS2_OUT_BASE"/*; do
    folder_name=$(basename "$folder")
    bed_file=$(find "$folder" -name "*_summits.bed" -type f | head -n 1)

    if [[ -f "$bed_file" ]]; then
        echo -e "${folder_name}\t${bed_file}" >> "$PEAK_LIST_FILE"
    else
        echo "⚠️  No narrowPeak file found for $folder_name"
    fi
done

# Run peak merging with R
echo "🔗 Merging peaks with iterative_overlap_peak_merging.R..."
Rscript "$MERGING_SCRIPT" \
    -i "$PEAK_LIST_FILE" \
    -g "$GENOME" \
    -d "$MERGED_PEAK_DIR"/ \
    -o "celltype_condition" \
    --blacklist "$BLACKLIST" \
    --chromSize "$CHROMSIZE"

echo "🧾 Logging R session info after merging..."
#Rscript -e "sink('results/logs/r_sessioninfo_peak_merging.txt'); sessionInfo(); sink()"

echo "✅ Peak calling and merging complete."

# Clean and sort union peakset
tail -n +2 "$UNION_PEAKSET" | awk '{OFS="\t"; print $1, $2, $3}' | sort -k1,1 -k2,2n > "$UNION_SORTED"
awk '{OFS="\t"; print $1, $2, $3, $1"-"$2"-"$3}' "$UNION_SORTED" > "results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted_label.peakset"

UNION_SORTED_LABEL_NEW="results/Human_cells/D3_ATAC/D3_merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted_label.peakset"
UNION_SORTED_LABEL_OLD="results/Human_cells/D3_ATAC/mergedpeaks/condition_treatment.filteredNfixed.union_clean_sorted_label.peakset"

OLD="$UNION_SORTED_LABEL_OLD"
NEW="$UNION_SORTED_LABEL_NEW"
MASTER_V2="results/Human_cells/D3_ATAC/D3_merged_peaks/master_v2.union.peakSet"

tmp_old=$(mktemp); tmp_new=$(mktemp)

# drop header lines, then exact-coordinate dedup across both sets
tail -n +2 "$OLD" > "$tmp_old"
tail -n +2 "$NEW" > "$tmp_new"

cat "$tmp_old" "$tmp_new" \
  | sort -k1,1 -k2,2n -k3,3n -u \
  > "$MASTER_V2"

rm -f "$tmp_old" "$tmp_new"

# sanity
wc -l "$OLD" "$NEW" "$MASTER_V2"
head -n 3 "$MASTER_V2


