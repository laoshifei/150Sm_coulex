#!/usr/bin/env bash
# Schedule native C++ extraction with bounded concurrency and complete waits.
set -euo pipefail
if (( $# < 1 || $# > 5 )); then
  echo "Usage: $0 OUTPUT_DIR [RUN_LIST] [JOBS=2] [SAMPLE_ENTRIES=2000000] [INPUT_DIR=/mnt/e/Sm150/rootfile]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
project="$(dirname -- "$script_dir")"
output="$1"
run_list="${2:-$project/pars/run.txt}"
jobs="${3:-2}"
sample="${4:-2000000}"
input="${5:-/mnt/e/Sm150/rootfile}"
[[ "$jobs" =~ ^[1-8]$ && "$sample" =~ ^[0-9]+$ ]] || { echo "Invalid jobs/sample count" >&2; exit 2; }
[[ -x "$script_dir/extract" ]] || { echo "Run make first" >&2; exit 2; }
mkdir -p -- "$output/logs"
pids=()
failed=0
while IFS= read -r run || [[ -n "$run" ]]; do
  [[ -z "$run" ]] && continue
  [[ "$run" =~ ^[0-9]{1,4}$ ]] || { echo "Invalid run: $run" >&2; exit 2; }
  printf -v rn '%04d' "$((10#$run))"
  "$script_dir/extract" "$((10#$run))" --project "$project" \
    --input "$input/data_C2_${rn}_${rn}.root" --output "$output/spectra_${rn}.root" \
    --sample-entries "$sample" > "$output/logs/extract_${rn}.log" 2>&1 &
  pids+=("$!")
  if (( ${#pids[@]} >= jobs )); then
    if ! wait "${pids[0]}"; then failed=1; fi
    pids=("${pids[@]:1}")
  fi
done < "$run_list"
for pid in "${pids[@]}"; do if ! wait "$pid"; then failed=1; fi; done
exit "$failed"
