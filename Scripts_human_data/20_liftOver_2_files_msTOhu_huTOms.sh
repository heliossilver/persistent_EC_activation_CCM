#!/usr/bin/env bash
set -euo pipefail

# ------------------------------
# Liftover peaks between genomes with reciprocal check
# Supports hg38 <-> mm10
# ------------------------------

usage() {
  cat <<'EOF'
Usage:
  liftover_peaks.sh \
    -a <PEAKS_A.bed[.gz]> -A <hg38|mm10> \
    -b <PEAKS_B.bed[.gz]> -B <hg38|mm10> \
    -r <reference_dir_with_chain_files> \
    -o <output_dir> \
    [--chainA2B <chainfile>] [--chainB2A <chainfile>] \
    [--minMatch 0.5] [--maxLen 1000] [--prefix NAME] [-v]

Notes:
- Requires: liftOver, bedtools, awk, zcat/gzip
- If the 4th BED column (name) is missing, it will be added as "chr-start-end"
- Filters to standard chromosomes and by length (<= maxLen) after liftover
- Produces reciprocal-filtered lifted sets and intersects them with the opposite species peakset

Outputs (in -o):
  NAME_AtoB.filtered.bed
  NAME_AtoB.reciprocal.bed
  NAME_BtoA.filtered.bed
  NAME_BtoA.reciprocal.bed
  NAME_AtoB_vs_B.intersect.tsv
  NAME_BtoA_vs_A.intersect.tsv
  unmapped_* files for diagnostics
EOF
  exit 1
}

# -------- args --------
PEAKS_A=""
PEAKS_B=""
GENOME_A=""
GENOME_B=""
REFDIR="reference"
OUTDIR="results/human_mouse_integration/liftover"
CHAIN_A2B=""
CHAIN_B2A=""
MINMATCH="0.5"
MAXLEN="1000"
PREFIX="liftover"
VERBOSE=0

vprint(){ [[ $VERBOSE -eq 1 ]] && echo "$@" || true; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) PEAKS_A="$2"; shift 2 ;;
    -b) PEAKS_B="$2"; shift 2 ;;
    -A) GENOME_A="$2"; shift 2 ;;
    -B) GENOME_B="$2"; shift 2 ;;
    -r) REFDIR="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    --chainA2B) CHAIN_A2B="$2"; shift 2 ;;
    --chainB2A) CHAIN_B2A="$2"; shift 2 ;;
    --minMatch) MINMATCH="$2"; shift 2 ;;
    --maxLen)   MAXLEN="$2"; shift 2 ;;
    --prefix)   PREFIX="$2"; shift 2 ;;
    -v) VERBOSE=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown arg: $1"; usage ;;
  esac
done

# -------- checks --------
[[ -n "$PEAKS_A" && -n "$PEAKS_B" && -n "$GENOME_A" && -n "$GENOME_B" ]] || usage
command -v liftOver >/dev/null 2>&1 || { echo "ERROR: liftOver not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
[[ -d "$REFDIR" ]] || { echo "ERROR: ref dir not found: $REFDIR"; exit 1; }

# auto-pick chains if not provided
auto_chain() {
  local src="$1" dst="$2"
  if [[ "$src" == "hg38" && "$dst" == "mm10" ]]; then
    echo "$REFDIR/hg38ToMm10.over.chain"
  elif [[ "$src" == "mm10" && "$dst" == "hg38" ]]; then
    echo "$REFDIR/mm10ToHg38.over.chain"
  else
    echo ""
    return 1
  fi
}

if [[ -z "$CHAIN_A2B" ]]; then
  CHAIN_A2B=$(auto_chain "$GENOME_A" "$GENOME_B") || { echo "ERROR: no chain for $GENOME_A->$GENOME_B"; exit 1; }
fi
if [[ -z "$CHAIN_B2A" ]]; then
  CHAIN_B2A=$(auto_chain "$GENOME_B" "$GENOME_A") || { echo "ERROR: no chain for $GENOME_B->$GENOME_A"; exit 1; }
fi
[[ -f "$CHAIN_A2B" ]] || { echo "ERROR: chain not found: $CHAIN_A2B"; exit 1; }
[[ -f "$CHAIN_B2A" ]] || { echo "ERROR: chain not found: $CHAIN_B2A"; exit 1; }

mkdir -p "$OUTDIR"

# std-chr regex per genome
stdchr_regex() {
  case "$1" in
    hg38) echo '^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)$' ;;   # chr1-22,X,Y
    mm10) echo '^chr([1-9]$|1[0-9]$|X$|Y$)$' ;;           # chr1-19,X,Y
    *) echo '.*' ;;
  esac
}

REGEX_A=$(stdchr_regex "$GENOME_A")
REGEX_B=$(stdchr_regex "$GENOME_B")

# -------------- helpers --------------
# read BED(.gz) and ensure >=3 cols, add 4th name if missing
normalize_bed() {
  local infile="$1" out="$2"
  vprint "Normalizing: $infile -> $out"
  if [[ "$infile" =~ \.gz$ ]]; then
    zcat -- "$infile"
  else
    cat -- "$infile"
  fi | awk -v OFS="\t" '
    BEGIN{ }
    $0 ~ /^#/ { next }
    NF>=3 {
      chr=$1; start=$2; end=$3;
      if (NF>=4 && $4!="") name=$4; else name=chr"-"start"-"end;
      # rebuild line with at least 4 fields; keep extras if present
      $1=chr; $2=start; $3=end; $4=name;
      print $0
    }' > "$out"
}

filter_len_chr() {
  local infile="$1" maxlen="$2" regex="$3" out="$4"
  awk -v OFS="\t" -v MAXLEN="$maxlen" -v RGX="$regex" '
    ($3-$2)<=MAXLEN && $1 ~ RGX' "$infile" > "$out"
}

# -------------- workflow --------------
TMP_A="$OUTDIR/${PREFIX}_A.norm.bed"
TMP_B="$OUTDIR/${PREFIX}_B.norm.bed"
normalize_bed "$PEAKS_A" "$TMP_A"
normalize_bed "$PEAKS_B" "$TMP_B"

# A->B
echo "▶️  Liftover A($GENOME_A) → B($GENOME_B)..."
liftOver "$TMP_A" "$CHAIN_A2B" \
  "$OUTDIR/${PREFIX}_AtoB.bed" \
  "$OUTDIR/${PREFIX}_unmapped_AtoB.bed" \
  -minMatch="$MINMATCH"

echo "🔎  Filter A->B by len≤$MAXLEN & std chroms of $GENOME_B..."
filter_len_chr "$OUTDIR/${PREFIX}_AtoB.bed" "$MAXLEN" "$REGEX_B" "$OUTDIR/${PREFIX}_AtoB.filtered.bed"

# B->A
echo "▶️  Liftover B($GENOME_B) → A($GENOME_A)..."
liftOver "$TMP_B" "$CHAIN_B2A" \
  "$OUTDIR/${PREFIX}_BtoA.bed" \
  "$OUTDIR/${PREFIX}_unmapped_BtoA.bed" \
  -minMatch="$MINMATCH"

echo "🔎  Filter B->A by len≤$MAXLEN & std chroms of $GENOME_A..."
filter_len_chr "$OUTDIR/${PREFIX}_BtoA.bed" "$MAXLEN" "$REGEX_A" "$OUTDIR/${PREFIX}_BtoA.filtered.bed"

# Reciprocal: lift back and keep names that return
echo "🔁  Reciprocal checks..."
# A->B back to A
liftOver "$OUTDIR/${PREFIX}_AtoB.filtered.bed" "$CHAIN_B2A" \
  "$OUTDIR/${PREFIX}_AtoB.back2A.bed" \
  "$OUTDIR/${PREFIX}_unmapped_AtoB_back.bed" \
  -minMatch="$MINMATCH"

# B->A back to B
liftOver "$OUTDIR/${PREFIX}_BtoA.filtered.bed" "$CHAIN_A2B" \
  "$OUTDIR/${PREFIX}_BtoA.back2B.bed" \
  "$OUTDIR/${PREFIX}_unmapped_BtoA_back.bed" \
  -minMatch="$MINMATCH"

# Keep only entries whose $4 name came back (ID match) and on std chroms
awk -v RGX="$REGEX_B" 'NR==FNR{ ids[$4]; next } ($4 in ids) && ($1 ~ RGX)' \
  "$OUTDIR/${PREFIX}_AtoB.back2A.bed" \
  "$OUTDIR/${PREFIX}_AtoB.filtered.bed" > "$OUTDIR/${PREFIX}_AtoB.reciprocal.bed"

awk -v RGX="$REGEX_A" 'NR==FNR{ ids[$4]; next } ($4 in ids) && ($1 ~ RGX)' \
  "$OUTDIR/${PREFIX}_BtoA.back2B.bed" \
  "$OUTDIR/${PREFIX}_BtoA.filtered.bed" > "$OUTDIR/${PREFIX}_BtoA.reciprocal.bed"

# Intersections with the opposite species’ original peaks
echo "🔗  Intersections with opposite peakset..."
bedtools intersect -wo \
  -a "$OUTDIR/${PREFIX}_AtoB.reciprocal.bed" \
  -b "$TMP_B" > "$OUTDIR/${PREFIX}_AtoB_vs_B.intersect.tsv"

bedtools intersect -wo \
  -a "$OUTDIR/${PREFIX}_BtoA.reciprocal.bed" \
  -b "$TMP_A" > "$OUTDIR/${PREFIX}_BtoA_vs_A.intersect.tsv"

echo "✅ Done. Outputs in: $OUTDIR"
