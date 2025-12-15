#!/bin/bash
set -e  # stop if any command fails

# === CONFIGURATION ===
BED_INPUT_DIR="results/ATAC/bedfiles_signac"
BLACKLIST="reference/mm10-blacklist.v2.bed.gz"
CHROMSIZE="reference/mm10.chrom.sizes"
BAM_FILE_DIR="results/ATAC/bam_files/"
BW_FILE_DIR="results/ATAC/bigwigs/"
MACS2_OUT_BASE="results/ATAC/macs2_peaks"
MERGED_PEAK_DIR="results/ATAC/merged_peaks/"
PEAK_LIST_FILE="results/ATAC/folder_summit_paths.txt"
MERGING_SCRIPT="scripts/utils/iterative_overlap_peak_merging.R"
GENOME="mm10"
CONDA_ENV_NAME="macs2_env"

# Create output folders
mkdir -p "$BAM_FILE_DIR"
mkdir -p "$BW_FILE_DIR"
mkdir -p "$MACS2_OUT_BASE"
mkdir -p "$MERGED_PEAK_DIR"
mkdir -p "$(dirname "$PEAK_LIST_FILE")"
mkdir -p "results/logs"

# === Activate conda ===
echo "⚙️  Activating conda environment: $CONDA_ENV_NAME"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV_NAME"

# Log MACS2 env info
echo "📝 Logging MACS2 conda environment info..."
conda info --envs > logs/macs2_env_log.txt
conda list > logs/macs2_env_packages.txt


# Iterate over each BED file in the peak_data_bed directory
for bed_file in "$BED_INPUT_DIR"/*.bed; do
    # Extract the filename without extension
    filename=$(basename "$bed_file" .bed)
    
    # Convert BED to BAM
    bedToBam -i "$bed_file" -g "$CHROMSIZE" > "$BAM_FILE_DIR""$filename".bam
    
    # Sort BAM
    samtools sort -o "$BAM_FILE_DIR""$filename"_sorted.bam "$BAM_FILE_DIR""$filename".bam
    
    # Index sorted BAM
    samtools index "$BAM_FILE_DIR""$filename"_sorted.bam
    
    # Calculate coverage and generate bigWig
    bamCoverage --normalizeUsing RPKM -b "$BAM_FILE_DIR""$filename"_sorted.bam -o "$BW_FILE_DIR""$filename".bw -p max
done


echo "🔍 Running MACS2 peak calling on each BED file..."

# Run MACS2 on each BED file
for file in "$BED_INPUT_DIR"/*.bed; do
    type=$(basename "$file" .bed)
    out_dir="$MACS2_OUT_BASE/$type"
    mkdir -p "$out_dir"

    echo "📈 Calling peaks for $type"
    macs2 callpeak \
        --treatment "$file" \
        --shift -75 \
        --extsize 150 \
        -g mm \
        --outdir "$out_dir" \
        --name "$type" \
        -q 0.05 \
        --call-summits \
        --nomodel \
        -f BED
done

# (Optional but recommended) deactivate conda before switching to R
echo "🚪 Deactivating macs2_env before iteration..."
conda deactivate


# Generate list of .narrowPeak files
echo "🧾 Creating peak list file: $PEAK_LIST_FILE"
echo -n "" > "$PEAK_LIST_FILE"


for folder in "$MACS2_OUT_BASE"/*; do
    folder_name=$(basename "$folder")
    bed_file=$(find "$folder" -name "*.bed" -type f | head -n 1)

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

