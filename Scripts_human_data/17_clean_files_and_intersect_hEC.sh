#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#                     SETTINGS
# ============================================================

# Overlap stringency toggle
F_ARGS=()                     # ≥1 bp overlap
# F_ARGS=( -f 0.5 -r )        # ≥50% reciprocal overlap (recommended for peaks)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

stdsort () {
  awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/' "$1" | sort -k1,1 -k2,2n
}

require() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found"; exit 1; }
}
require bedtools

# ============================================================
#              MAKE CONSENSUS FROM REPLICATES
# ============================================================

make_consensus () {
  local BASE_DIR="$1"
  local UNION="$2"
  local OUTFILE="$3"
  shift 3
  local REPS=("$@")

  # start from union
  local TMP="${OUTFILE}.tmp"
  cp "$UNION" "$TMP"

  for rep in "${REPS[@]}"; do
    bedtools intersect -a "$TMP" -b "$rep" -u "${F_ARGS[@]}" > "${TMP}.next"
    mv "${TMP}.next" "$TMP"
  done

  mv "$TMP" "$OUTFILE"
  echo "   - wrote $OUTFILE"
}

# ============================================================
#                     PROCESS HPAEC
# ============================================================

process_hpaec () {
  BASE_DIR="results/Human_cells/HPAEC_ATAC/HPAEC_merged_peaks"
  UNION="${BASE_DIR}/union.std.sorted.bed"

  # prep inputs
  stdsort "${BASE_DIR}/celltype_condition.filteredNfixed.union_clean_sorted.peakset" > "$UNION"

  for r in V1 V2 C1 C2; do
    stdsort "${BASE_DIR}/HPAEC_${r}.filterNfixed.peakset" > "${BASE_DIR}/${r}.std.sorted.bed"
  done

  make_consensus "$BASE_DIR" "$UNION" "${BASE_DIR}/vehicle_consensus.bed" \
    "${BASE_DIR}/V1.std.sorted.bed" "${BASE_DIR}/V2.std.sorted.bed"

  make_consensus "$BASE_DIR" "$UNION" "${BASE_DIR}/ccm_consensus.bed" \
    "${BASE_DIR}/C1.std.sorted.bed" "${BASE_DIR}/C2.std.sorted.bed"
}

# ============================================================
#                     PROCESS HUVEC
# ============================================================

process_huvec () {
  BASE_DIR="results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks"
  UNION="${BASE_DIR}/union.std.sorted.bed"

  stdsort "${BASE_DIR}/celltype_condition.filteredNfixed.union_clean_sorted.peakset" > "$UNION"

  for r in V1 V2 C1 C2; do
    stdsort "${BASE_DIR}/HUVEC_${r}.filterNfixed.peakset" > "${BASE_DIR}/${r}.std.sorted.bed"
  done

  make_consensus "$BASE_DIR" "$UNION" "${BASE_DIR}/vehicle_consensus.bed" \
    "${BASE_DIR}/V1.std.sorted.bed" "${BASE_DIR}/V2.std.sorted.bed"

  make_consensus "$BASE_DIR" "$UNION" "${BASE_DIR}/ccm_consensus.bed" \
    "${BASE_DIR}/C1.std.sorted.bed" "${BASE_DIR}/C2.std.sorted.bed"
}

# ============================================================
#                     PROCESS D3 (triplicates)
# ============================================================

process_d3 () {
  BASE_DIR="results/Human_cells/D3_ATAC/D3_merged_peaks"
  UNION="${BASE_DIR}/union.std.sorted.bed"

  # filter union once
  stdsort "${BASE_DIR}/celltype_condition.filteredNfixed.union_clean_sorted.peakset" > "$UNION"

  # filter & sort replicates
  for r in WT_V1 WT_V2 WT_V3 WT_C1 WT_C2 WT_C3; do
    stdsort "${BASE_DIR}/${r}.filterNfixed.peakset" > "${BASE_DIR}/${r}.std.sorted.bed"
  done

  # Vehicle consensus (all 3 must overlap)
  make_consensus "$BASE_DIR" "$UNION" "${BASE_DIR}/vehicle_consensus.bed" \
    "${BASE_DIR}/WT_V1.std.sorted.bed" \
    "${BASE_DIR}/WT_V2.std.sorted.bed" \
    "${BASE_DIR}/WT_V3.std.sorted.bed"

  # Control consensus (all 3 must overlap)
  make_consensus "$BASE_DIR" "$UNION" "${BASE_DIR}/ccm_consensus.bed" \
    "${BASE_DIR}/WT_C1.std.sorted.bed" \
    "${BASE_DIR}/WT_C2.std.sorted.bed" \
    "${BASE_DIR}/WT_C3.std.sorted.bed"
}

# ============================================================
#                     RUN ALL
# ============================================================

process_hpaec
process_huvec
process_d3

# ============================================================
#              VENN HELPER FOR 3 SETS (A, B, C)
# ============================================================

venn_for_three () {
  local A_STD="$1"    # path to set A bed
  local B_STD="$2"    # path to set B bed
  local C_STD="$3"    # path to set C bed
  local OUT_DIR="$4"  # output dir
  local PREFIX="$5"   # label (e.g., vehicle or ccm)

  mkdir -p "$OUT_DIR"

  # ---- count directional overlaps ----
  dir_count()  { bedtools intersect -u -a "$1" -b "$2" "${F_ARGS[@]}" | wc -l; }
  triple_from(){ bedtools intersect -u -a "$1" -b "$2" "${F_ARGS[@]}" \
                 | bedtools intersect -u -a - -b "$3" "${F_ARGS[@]}" | wc -l; }

  # ---- areas ----
  area1=$(wc -l < "$A_STD")
  area2=$(wc -l < "$B_STD")
  area3=$(wc -l < "$C_STD")

  # ---- pairwise (symmetric via min of directions) ----
  n12_A=$(dir_count "$A_STD" "$B_STD")
  n12_B=$(dir_count "$B_STD" "$A_STD")
  n12=$(( n12_A < n12_B ? n12_A : n12_B ))

  n13_A=$(dir_count "$A_STD" "$C_STD")
  n13_C=$(dir_count "$C_STD" "$A_STD")
  n13=$(( n13_A < n13_C ? n13_A : n13_C ))

  n23_B=$(dir_count "$B_STD" "$C_STD")
  n23_C=$(dir_count "$C_STD" "$B_STD")
  n23=$(( n23_B < n23_C ? n23_B : n23_C ))

  # ---- triple (symmetric via min of 3 perspectives) ----
  n123_A=$(triple_from "$A_STD" "$B_STD" "$C_STD")
  n123_B=$(triple_from "$B_STD" "$A_STD" "$C_STD")
  n123_C=$(triple_from "$C_STD" "$A_STD" "$B_STD")
  n123=$(( n123_A < n123_B ? n123_A : n123_B ))
  n123=$(( n123 < n123_C ? n123 : n123_C ))

  # ---- write counts ----
  VENN_OUT="${OUT_DIR}/venn_counts_${PREFIX}.tsv"
  {
    printf "area1\t%s\n" "$area1"
    printf "area2\t%s\n" "$area2"
    printf "area3\t%s\n" "$area3"
    printf "n12\t%s\n"   "$n12"
    printf "n13\t%s\n"   "$n13"
    printf "n23\t%s\n"   "$n23"
    printf "n123\t%s\n"  "$n123"
  } > "$VENN_OUT"

  echo "Wrote: $VENN_OUT"

  # ---- optional BED exports ----
  bedtools intersect -u -a "$A_STD" -b "$B_STD" "${F_ARGS[@]}" \
  | bedtools intersect -u -a - -b "$C_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/triple_${PREFIX}_fromA.bed"

  bedtools intersect -u -a "$A_STD" -b "$B_STD" "${F_ARGS[@]}" \
  | bedtools intersect -v -a - -b "$C_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/AB_only_${PREFIX}_fromA.bed"

  bedtools intersect -u -a "$A_STD" -b "$C_STD" "${F_ARGS[@]}" \
  | bedtools intersect -v -a - -b "$B_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/AC_only_${PREFIX}_fromA.bed"

  bedtools intersect -u -a "$B_STD" -b "$C_STD" "${F_ARGS[@]}" \
  | bedtools intersect -v -a - -b "$A_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/BC_only_${PREFIX}_fromB.bed"

  bedtools intersect -v -a "$A_STD" -b "$B_STD" "${F_ARGS[@]}" \
  | bedtools intersect -v -a - -b "$C_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/A_unique_${PREFIX}.bed"

  bedtools intersect -v -a "$B_STD" -b "$A_STD" "${F_ARGS[@]}" \
  | bedtools intersect -v -a - -b "$C_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/B_unique_${PREFIX}.bed"

  bedtools intersect -v -a "$C_STD" -b "$A_STD" "${F_ARGS[@]}" \
  | bedtools intersect -v -a - -b "$B_STD" "${F_ARGS[@]}" \
  > "${OUT_DIR}/C_unique_${PREFIX}.bed"
}

# ============================================================
#              RUN VENN FOR VEHICLE AND CCM
# ============================================================

VENN_OUTDIR="results/Human_cells/common/venn"

# VEHICLE: HPAEC, HUVEC, D3
venn_for_three \
  "results/Human_cells/HPAEC_ATAC/HPAEC_merged_peaks/vehicle_consensus.bed" \
  "results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks/vehicle_consensus.bed" \
  "results/Human_cells/D3_ATAC/D3_merged_peaks/vehicle_consensus.bed" \
  "$VENN_OUTDIR" "vehicle"

# CCM: HPAEC, HUVEC, D3
venn_for_three \
  "results/Human_cells/HPAEC_ATAC/HPAEC_merged_peaks/ccm_consensus.bed" \
  "results/Human_cells/HUVEC_ATAC/HUVEC_merged_peaks/ccm_consensus.bed" \
  "results/Human_cells/D3_ATAC/D3_merged_peaks/ccm_consensus.bed" \
  "$VENN_OUTDIR" "ccm"
