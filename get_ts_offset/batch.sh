#!/usr/bin/env bash
# Specify the run list, input, output and concurrency explicitly; wait for every job and collect exit codes.
set -euo pipefail
if (( $# < 3 || $# > 5 )); then
  echo 'Usage: batch.sh RUN_LIST INPUT_DIR OUTPUT_DIR [JOBS=2] [raw|residual]' >&2
  exit 2
fi
list=$1
input_dir=$(realpath -- "$2")
output_dir=$(realpath -m -- "$3")
jobs_limit=${4:-2}
mode=${5:-raw}
[[ $jobs_limit =~ ^[1-9][0-9]*$ && $mode =~ ^(raw|residual)$ ]] || exit 2
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
[[ -x "$script_dir/analysis" ]] || { echo 'Run make -C get_ts_offset first' >&2; exit 2; }
[[ -f "$list" ]] || { echo "Missing run list: $list" >&2; exit 2; }
mkdir -p -- "$(dirname -- "$output_dir")"
mkdir -- "$output_dir"
declare -A active=() seen=()
failed=0
reap() {
  local pid rc=0
  wait -n -p pid "${!active[@]}" || rc=$?
  printf '%s\t%s\n' "${active[$pid]}" "$rc" >> "$output_dir/status.tsv"
  unset 'active[$pid]'
  if ((rc)); then failed=1; fi
}
while IFS= read -r line || [[ -n $line ]]; do
  line=${line%%#*}
  read -r run extra <<< "$line"
  [[ -n ${run:-} ]] || continue
  [[ $run =~ ^[0-9]{1,4}$ && -z ${extra:-} ]] || { echo "Invalid run row: $line" >&2; failed=1; continue; }
  printf -v run '%04d' "$((10#$run))"
  [[ -z ${seen[$run]+x} ]] || { echo "Duplicate run: $run" >&2; failed=1; continue; }
  seen[$run]=1
  "$script_dir/analysis" "$run" --mode "$mode" \
    --input "$input_dir/data_C2_${run}_${run}.root" \
    --output "$output_dir/data${run}_timing_${mode}.root" \
    > "$output_dir/${run}.log" 2>&1 &
  active[$!]=$run
  while ((${#active[@]} >= jobs_limit)); do reap; done
done < "$list"
while ((${#active[@]})); do reap; done
exit "$failed"
