#!/usr/bin/env bash
# Explicit inputs and bounded concurrency; no implicit full-project queue.
set -euo pipefail
if [[ $# -ne 5 ]]; then
  echo "Usage: $0 RUN_LIST INPUT_DIR NEW_OUTPUT_DIR JOBS SAMPLE_ENTRIES (0 = full input)" >&2
  exit 2
fi
run_list=$(realpath "$1")
input_dir=$(realpath "$2")
output_dir=$(realpath -m "$3")
jobs=$4
sample=$5
project=$(cd "$(dirname "$0")/.." && pwd)
[[ $jobs =~ ^[1-9][0-9]*$ && $jobs -le 64 && $sample =~ ^[0-9]+$ ]] || exit 2
[[ ! -e $output_dir && ! -L $output_dir ]] || { echo "Output directory already exists" >&2; exit 2; }
mapfile -t runs < <(awk 'NF && $1 !~ /^#/ {print $1}' "$run_list")
[[ ${#runs[@]} -gt 0 ]] || { echo "Empty run list" >&2; exit 2; }
declare -A seen=()
for run in "${runs[@]}"; do
  [[ $run =~ ^[0-9]{1,4}$ ]] || { echo "Invalid run: $run" >&2; exit 2; }
  printf -v key '%04d' "$((10#$run))"
  [[ -z ${seen[$key]:-} ]] || { echo "Duplicate run: $key" >&2; exit 2; }
  seen[$key]=1
  [[ -f "$input_dir/s3_$key.root" ]] || { echo "Missing s3_events input: $key" >&2; exit 2; }
done
mkdir -p "$output_dir/logs"
active=0
failed=0
for run in "${runs[@]}"; do
  printf -v key '%04d' "$((10#$run))"
  "$project/s3_coin_ana/ana" "$key" --project "$project" \
    --input "$input_dir/s3_$key.root" --output "$output_dir/s3_$key.root" \
    --sample-entries "$sample" --segments 20 > "$output_dir/logs/$key.log" 2>&1 &
  active=$((active+1))
  if (( active >= jobs )); then
    if ! wait -n; then active=$((active-1)); failed=1; break; fi
    active=$((active-1))
  fi
done
while (( active > 0 )); do
  if ! wait -n; then failed=1; fi
  active=$((active-1))
done
exit "$failed"
