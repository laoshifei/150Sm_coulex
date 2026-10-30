#!/usr/bin/env bash
set -euo pipefail
if (( $# != 2 )); then
  echo "Usage: $0 SPECTRA_ROOT NEW_FIGURE_DIR" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
# Reject characters that would alter the ROOT expression.
for arg in "$1" "$2"; do
  [[ "$arg" != *'"'* && "$arg" != *'\'* && "$arg" != *$'\n'* ]] || { echo "Unsupported path characters" >&2; exit 2; }
done
exec root -l -b -q -n "$script_dir/draw_sample.cpp(\"$1\",\"$2\")"
