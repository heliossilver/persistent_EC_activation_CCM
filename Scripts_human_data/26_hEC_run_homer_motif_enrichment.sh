#!/usr/bin/env bash
set -euo pipefail

# -------- helpers --------
cpu_count() {
  if command -v nproc >/dev/null 2>&1; then nproc
  elif [[ "$(uname -s)" == "Darwin" ]]; then sysctl -n hw.ncpu
  else echo 1
  fi
}

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found"; exit 1; }; }

run_enrichment() {
  local input_dir="$1"
  local output_dir="$2"
  local bg_file="$3"
  local genome="${4:-hg38}"     # optional 4th arg (defaults to hg38)
  local size="${5:-200}"        # optional 5th arg (defaults to 200)
  local de_novo="${6:-no}"      # "yes" to allow de-novo, "no" to pass -nomotif

  echo "🚀 Running enrichment in: $input_dir"
  echo "   ➤ Output to:          $output_dir"
  echo "   ➤ Background:         $bg_file"
  echo "   ➤ Genome/size:        $genome / $size"
  echo "   ➤ De-novo motifs:     $de_novo"

  [[ -d "$input_dir" ]] || { echo "ERROR: input_dir not found: $input_dir"; exit 1; }
  [[ -f "$bg_file"   ]] || { echo "ERROR: background file not found: $bg_file"; exit 1; }
  mkdir -p "$output_dir"

  # build HOMER de-novo flag
  local nomotif_flag=()
  [[ "$de_novo" == "no" ]] && nomotif_flag=(-nomotif)

  # Find .bed files robustly (handles spaces) and run in parallel
  find "$input_dir" -maxdepth 1 -type f -name "*.bed" -print0 |
  parallel -0 --jobs "$(cpu_count)" --bar '
    file_path="{}"
    filename="$(basename "$file_path" .bed)"
    outdir="'"$output_dir"'/${filename}"
    mkdir -p "$outdir"
    echo "🔍 Processing: $file_path → $outdir"
    findMotifsGenome.pl "$file_path" "'"$genome"'" "$outdir" -size "'"$size"'" -bg "'"$bg_file"'" '"${nomotif_flag[@]}"'
    echo "✅ Completed:  $filename"
  '
}

# -------- preflight --------
need find
need parallel
need findMotifsGenome.pl

# -------- runs --------
# 1) hEC response 
run_enrichment \
  "results/Human_cells/Homer/hEC_response_CCM" \
  "results/Human_cells/Homer/hEC_response_CCM_enrichment" \
  "results/Human_cells/data/master_univers.peakset" \
  "hg38" 200 "no"

echo "🎉 All motif enrichment analyses are complete!"

