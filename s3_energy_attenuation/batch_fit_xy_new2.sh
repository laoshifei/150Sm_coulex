#!/usr/bin/env bash
set -euo pipefail
if (( $# < 2 || $# > 3 )); then
  echo "Usage: $0 SPECTRA_DIR NEW_RESULT_DIR [RUN_LIST]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
project="$(dirname -- "$script_dir")"
exec "$script_dir/calibrate" --input-dir "$1" --output-dir "$2" --run-list "${3:-$project/pars/run.txt}"

