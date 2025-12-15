#!/bin/bash
set -euo pipefail  # safer shell

# === CONFIGURATION ===
BLACKLIST="reference/hg38-blacklist.v2.bed.gz"
CHROMSIZE="reference/hg38.chrom.sizes"
OUT_BASE="results/Human_cells/HUVEC_ATAC"
MERGED_PEAK_DIR="${OUT_BASE}/HUVEC_merged_peaks"
PEAK_LIST_FILE="${OUT_BASE}/folder_summit_paths.txt"
MACS2_OUT_BASE="${OUT_BASE}/macs2"
MERGING_SCRIPT="scripts/utils/iterative_overlap_peak_merging.R"
GENOME="hg38"

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

