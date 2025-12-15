#!/bin/bash

set -euo pipefail

# Change directory to where the peak files are located
cd results/Human_cells/D3_ATAC/D3_merged_peaks || { echo "cd failed"; exit 1; }


# Required union file (make sure it exists here)
UNION="celltype_condition.filteredNfixed.union_clean_sorted.peakset"
[[ -s "$UNION" ]] || { echo "ERROR: Missing union file: $UNION"; exit 1; }


# Create directories to store the combined peaksets and intersected peaksets
# Create directories to store the combined peaksets and intersected peaksets
mkdir -p combined_peak_by_condition
mkdir -p peakset_by_celltype
mkdir -p peakset_by_comparisons

# List of cell types
celltypes=("KD_C" "KD_V" "WT_C" "WT_V" "KD_T" "WT_T" "KC_T" "WC_T" "K1_C" "K1_V" "W1_C" "W1_V")

# Function to process cell type peaksets
# Function to process cell type peaksets
process_celltype() {
  celltype=$1
  echo "Processing $celltype..."

  # Combine all three replicates for the cell type
  > combined_peak_by_condition/${celltype}_combined.peakset  # Create or clear the file
  for replicate in 1 2 3; do
    # Check if the replicate file exists
    if [[ -f ${celltype}${replicate}.filterNfixed.peakset ]]; then
      tail -n +2 ${celltype}${replicate}.filterNfixed.peakset >> combined_peak_by_condition/${celltype}_combined.peakset
    else
      echo "Warning: File ${celltype}${replicate}.filterNfixed.peakset not found!"
    fi
  done

  # Extract the first three columns (seqnames, start, end) from the combined peakset
  awk '{OFS="\t"; print $1, $2, $3}' combined_peak_by_condition/${celltype}_combined.peakset > combined_peak_by_condition/${celltype}_combined_clean.peakset

  # Sort the cleaned combined peakset
  sort -k1,1 -k2,2n combined_peak_by_condition/${celltype}_combined_clean.peakset > combined_peak_by_condition/${celltype}_combined_clean_sorted.peakset

  # Perform bedtools intersect with the cleaned and sorted union peakset
  bedtools intersect -a combined_peak_by_condition/${celltype}_combined_clean_sorted.peakset -b celltype_condition.filteredNfixed.union_clean_sorted.peakset -wb > peakset_by_celltype/${celltype}.filteredNfixed.union_no_header.peakset

  # Extract only the last three columns (from the union peakset) and remove duplicates
  awk '{OFS="\t"; print $4, $5, $6}' peakset_by_celltype/${celltype}.filteredNfixed.union_no_header.peakset | sort -k1,1 -k2,2n | uniq > peakset_by_celltype/${celltype}.filteredNfixed.union_unique.peakset

  # Add the header back to the final file
  echo -e "chr\tstart\tend" | cat - peakset_by_celltype/${celltype}.filteredNfixed.union_unique.peakset > peakset_by_celltype/${celltype}.filteredNfixed.union.peakset

  # Clean up intermediate files
  rm -f peakset_by_celltype/${celltype}.filteredNfixed.union_no_header.peakset peakset_by_celltype/${celltype}.filteredNfixed.union_unique.peakset
  echo "$celltype processing completed!"
}

# Loop through each cell type
for celltype in "${celltypes[@]}"; do
  process_celltype "$celltype"
done

echo "All cell types processed!"


# List of comparisons
comparisons=("KD_V_vs_WT_V" "KD_V_vs_WT_C" "KD_C_vs_KD_V" "KD_C_vs_WT_V" "KD_C_vs_WT_C" "KD_T_vs_KD_V" "KD_T_vs_WT_V" "KD_T_vs_WT_T" "KC_T_vs_KD_C" "KC_T_vs_KD_T" "KC_T_vs_KD_V" "KC_T_vs_WT_V" "KC_T_vs_WT_T" "KC_T_vs_WC_T" "WT_C_vs_WT_V" "WT_T_vs_WT_V" "WC_T_vs_WT_V" "K1_C_vs_K1_V" "K1_C_vs_W1_V" "W1_C_vs_W1_V" "KD_C_vs_K1_C" "WT_C_vs_W1_C")

# Function to process comparisons

process_comparison() {
  comparison=$1
  echo "Processing comparison: $comparison..."

  # Extract cell types from the comparison name
  celltype1=$(echo "$comparison" | cut -d'_' -f1-2)
  celltype2=$(echo "$comparison" | cut -d'_' -f4-5)

  # Check if the source peakset files exist
  if [[ ! -f peakset_by_celltype/${celltype1}.filteredNfixed.union.peakset ]]; then
    echo "Error: File peakset_by_celltype/${celltype1}.filteredNfixed.union.peakset not found!"
    return
  fi

  if [[ ! -f peakset_by_celltype/${celltype2}.filteredNfixed.union.peakset ]]; then
    echo "Error: File peakset_by_celltype/${celltype2}.filteredNfixed.union.peakset not found!"
    return
  fi

  # Remove the header from both celltype peaksets
  tail -n +2 peakset_by_celltype/${celltype1}.filteredNfixed.union.peakset > tmp_1.peakset
  tail -n +2 peakset_by_celltype/${celltype2}.filteredNfixed.union.peakset > tmp_2.peakset

  # Check if tmp files were created successfully
  if [[ ! -s tmp_1.peakset ]]; then
    echo "Error: tmp_1.peakset is empty for $comparison!"
    return
  fi

  if [[ ! -s tmp_2.peakset ]]; then
    echo "Error: tmp_2.peakset is empty for $comparison!"
    return
  fi

  # Concatenate the two celltype peaksets
  cat tmp_1.peakset tmp_2.peakset > combined_peak_by_condition/${celltype1}_vs_${celltype2}_combined.peakset

  # Extract the first three columns
  awk '{OFS="\t"; print $1, $2, $3}' combined_peak_by_condition/${celltype1}_vs_${celltype2}_combined.peakset > combined_peak_by_condition/${celltype1}_vs_${celltype2}_combined_clean.peakset

  # Sort the cleaned peak set
  sort -k1,1 -k2,2n combined_peak_by_condition/${celltype1}_vs_${celltype2}_combined_clean.peakset > combined_peak_by_condition/${celltype1}_vs_${celltype2}_combined_clean_sorted.peakset

  # Perform bedtools intersect
  bedtools intersect -a combined_peak_by_condition/${celltype1}_vs_${celltype2}_combined_clean_sorted.peakset -b celltype_condition.filteredNfixed.union_clean_sorted.peakset -wb > peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_no_header.peakset

  # Check if the intersect file was created
  if [[ ! -s peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_no_header.peakset ]]; then
    echo "Warning: No intersected peaks found for comparison $comparison."
    return
  fi

  # Extract only the last three columns and remove duplicates
  awk '{OFS="\t"; print $4, $5, $6}' peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_no_header.peakset | sort -k1,1 -k2,2n | uniq > peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_unique.peakset

  # Add the header back if the file is not empty
  if [[ -s peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_unique.peakset ]]; then
    echo -e "chr\tstart\tend" | cat - peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_unique.peakset > peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union.peakset
  fi

  # Clean up intermediate files
  rm -f tmp_1.peakset tmp_2.peakset
  rm -f peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_no_header.peakset
  rm -f peakset_by_comparisons/${celltype1}_vs_${celltype2}.filteredNfixed.union_unique.peakset

  echo "$comparison processing completed!"
}

# Loop through each comparison
for comparison in "${comparisons[@]}"; do
  process_comparison "$comparison"
done


echo "All comparisons processed!"

