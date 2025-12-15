#!/bin/bash

# script12_run_homer_motif_enrichment.sh
# 🔍 This script runs HOMER motif enrichment for:
# - enhancer_flox
# - pathway_enhancers
# - enhancers_KO
# - cCRE_KO

set -e  # Stop if any command fails

# ===== Function: Enrichment Runner =====
run_enrichment() {
    input_dir=$1
    output_dir=$2
    bg_file=$3

    echo "🚀 Running enrichment in $input_dir using background $bg_file"

    mkdir -p "$output_dir"

    find "$input_dir" -maxdepth 1 -type f -name "*.txt" | parallel --jobs "$(nproc)" --bar '
        file_path={};
        filename=$(basename "$file_path" .txt);
        full_output_dir="'$output_dir'/${filename}";
        mkdir -p "$full_output_dir";
        tail -n +2 "$file_path" > "temp_${filename}.txt";
        echo "🔍 Processing $file_path -> $full_output_dir";
        findMotifsGenome.pl "temp_${filename}.txt" mm10 "$full_output_dir" -size 200 -bg "$bg_file" -nomotif;
        rm "temp_${filename}.txt";
        echo "✅ Completed $file_path"
    '
}

# ===== Run for enhancer_flox =====
run_enrichment "results/RNA_ATAC_integration/Homer/enhancer_flox" \
               "results/RNA_ATAC_integration/Homer/enhancer_flox_enrichment" \
               "results/ATAC/filtered_peaksets/peakset_all_flox/all_flox_filtered_union.peakset"

# ===== Run for pathway_enhancers =====
run_enrichment "results/RNA_ATAC_integration/Homer/pathway_enhancers" \
               "results/RNA_ATAC_integration/Homer/pathway_enhancers_enrichment" \
               "results/ATAC/filtered_peaksets/peakset_all_flox/all_flox_filtered_union.peakset"

# ===== Run for enhancers_KO =====
run_enrichment "results/RNA_ATAC_integration/Homer/enhancers_KO" \
               "results/RNA_ATAC_integration/Homer/enhancers_KO_enrichment" \
               "results/ATAC/merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted.peakset"

# ===== Run for cCRE_KO ===== 
run_enrichment "results/RNA_ATAC_integration/Homer/cCRE_KO" \
               "results/RNA_ATAC_integration/Homer/cCRE_KO_enrichment" \
               "results/ATAC/merged_peaks/celltype_condition.filteredNfixed.union_clean_sorted.peakset"

echo "🎉 All motif enrichment analyses are complete!"

