#!/usr/bin/env bash
# Run an explicit list with bounded concurrency; never replace existing results.
set -euo pipefail
if (( $# < 3 || $# > 4 )); then
  echo "Usage: $0 RUN_LIST INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY [JOBS=1]" >&2
  exit 2
fi
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(dirname -- "$script_dir")
run_list=$1
input_dir=$2
output_dir=$3
max_jobs=${4:-1}
[[ $max_jobs =~ ^[1-9][0-9]*$ ]] || { echo "Invalid job count" >&2; exit 2; }
[[ ! -e $output_dir ]] || { echo "Output directory already exists" >&2; exit 1; }
[[ -f $run_list && -d $input_dir ]] || { echo "Missing run list or input directory" >&2; exit 1; }
runs=()
declare -A seen=()
while IFS= read -r line || [[ -n $line ]]; do
  line=${line%%#*}
  read -r value extra <<< "$line"
  [[ -n $value ]] || continue
  [[ $value =~ ^[0-9]{1,4}$ && -z $extra ]] || { echo "Invalid run: $line" >&2; exit 2; }
  printf -v run '%04d' "$((10#$value))"
  [[ -z ${seen[$run]+x} ]] || { echo "Duplicate run: $run" >&2; exit 2; }
  seen[$run]=1
  [[ -r $input_dir/s3_${run}.root ]] || { echo "Missing input for $run" >&2; exit 1; }
  runs+=("$run")
done < "$run_list"
(( ${#runs[@]} > 0 )) || { echo "Empty run list" >&2; exit 2; }
mkdir -p -- "$output_dir/logs"
pids=()
status=0
for run in "${runs[@]}"; do
  "$script_dir/calibrate" --input "$input_dir/s3_${run}.root" --outdir "$output_dir/$run" > "$output_dir/logs/$run.log" 2>&1 &
  pids+=("$!")
  if (( ${#pids[@]} >= max_jobs )); then
    if ! wait "${pids[0]}"; then status=1; fi
    pids=("${pids[@]:1}")
    (( status == 0 )) || break
  fi
done
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
exit "$status"
