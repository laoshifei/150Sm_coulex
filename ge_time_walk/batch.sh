#!/usr/bin/env bash
# 原始读取在C++完成；显式名单/路径，持续补位并等待所有任务结束。
set -uo pipefail
if (( $# < 3 || $# > 5 )); then
  echo 'Usage: batch.sh RUN_LIST INPUT_DIR NEW_OUTPUT_DIR [JOBS=2] [PROFILE=local]' >&2
  exit 2
fi
list=$1
input=$2
output=$3
jobs_limit=${4:-2}
profile=${5:-local}
[[ $jobs_limit =~ ^[1-9][0-9]*$ ]] && (( jobs_limit <= 8 )) || { echo 'JOBS must be1..8' >&2; exit 2; }
directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project=$(cd -- "$directory/.." && pwd)
[[ -f $list && -d $input && -x $directory/build && ! -e $output ]] || { echo 'Check list/input/build; output directory must be new' >&2; exit 2; }
if [[ $profile == active ]]; then read -r profile < "$project/pars/time_alignment/active.txt"; fi
profile_args=(--profile "$profile")
if [[ $profile == local ]]; then
  # 整批固定同一目录；读取器仍在各轮核对完整参数SHA。
  profile_dir=$(python3 "$directory/parameters/manage.py" path) || exit 2
  profile_args=(--profile-dir "$profile_dir")
fi
declare -a runs=()
declare -A seen=() running=()
while IFS= read -r line || [[ -n $line ]]; do
  line=${line%%#*}
  [[ $line =~ ^[[:space:]]*$ ]] && continue
  [[ $line =~ ^[[:space:]]*([0-9]{1,4})[[:space:]]*$ ]] || { echo "Invalid run: $line" >&2; exit 2; }
  run=$((10#${BASH_REMATCH[1]}))
  (( run > 0 )) || exit 2
  printf -v label '%04d' "$run"
  [[ ! -v seen[$label] && -r $input/data_C2_${label}_${label}.root ]] || { echo "Duplicate/unreadable run: $label" >&2; exit 2; }
  seen[$label]=1;runs+=("$label")
done < "$list"
(( ${#runs[@]} > 0 )) || { echo 'Empty run list' >&2; exit 2; }
mkdir -- "$output" || exit 2
printf 'run\texit_code\n' > "$output/exit_status.tsv"
failures=0
finish_one() {
  local done_pid code
  wait -n -p done_pid "${!running[@]}" 2>/dev/null;code=$?
  # 已结束的子进程可能已离开jobs列表；指定PID的wait仍保存其退出状态。
  if [[ -z ${done_pid:-} ]]; then
    for done_pid in "${!running[@]}"; do wait "$done_pid";code=$?;break;done
  fi
  [[ -n ${done_pid:-} && -v running[$done_pid] ]] || { echo 'Cannot identify child status' >&2; exit 3; }
  printf '%s\t%s\n' "${running[$done_pid]}" "$code" >> "$output/exit_status.tsv"
  (( code == 0 )) || failures=$((failures+1))
  unset 'running[$done_pid]'
}
for label in "${runs[@]}"; do
  "$directory/build" "$((10#$label))" --project "$project" "${profile_args[@]}" \
    --input "$input/data_C2_${label}_${label}.root" --output "$output/$label.root" \
    > "$output/$label.log" 2>&1 &
  running[$!]=$label
  while (( ${#running[@]} >= jobs_limit )); do finish_one; done
done
while (( ${#running[@]} > 0 )); do finish_one; done
printf 'Finished %s runs; failures=%s\n' "${#runs[@]}" "$failures"
(( failures == 0 ))

