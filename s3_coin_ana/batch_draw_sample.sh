#!/usr/bin/env bash
# Render explicitly selected analysis files, one ROOT interpreter at a time.
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 NEW_FIGURE_DIR ANALYSIS.root [...]" >&2
  exit 2
fi
destination=$(realpath -m "$1")
shift
[[ ! -e $destination && ! -L $destination ]] || { echo "Output exists" >&2; exit 2; }
project=$(cd "$(dirname "$0")/.." && pwd)
for file in "$@"; do
  [[ -f $file ]] || { echo "Missing input: $file" >&2; exit 2; }
done
mkdir -p "$destination"
for file in "$@"; do
  input=$(realpath "$file")
  label=$(basename "$file" .root)
  # Quotes and backslashes are rejected before constructing a ROOT macro call.
  case "$input$destination$project" in *\"*|*\\*) echo "Unsupported path characters" >&2; exit 2;; esac
  root -l -b -q "$project/s3_coin_ana/draw_sample.cpp(\"$input\",\"$destination/$label\")"
done
