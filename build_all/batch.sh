#!/usr/bin/env bash
# ./batch.sh [并行数=8] [轮次名单=../pars/run.txt] [输出目录=项目rootfile]
# C++处理事件；本脚本只启动、补位和汇总退出状态。
set -uo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(dirname -- "$script_dir")"
max_jobs="${1:-8}"
run_list="${2:-$project_dir/pars/run.txt}"
output_dir="${3:-$project_dir/rootfile}"

if [[ $# -gt 3 || ! "$max_jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "用法：$0 [正整数并行数] [轮次名单] [输出目录]" >&2
  exit 2
fi
if [[ ! -x "$script_dir/build" || ! -r "$run_list" ]]; then
  echo "请先 make -C build_all，并检查轮次名单：$run_list" >&2
  exit 2
fi

# 核对整份名单和输出冲突后才启动任务。
runs=()
declare -A seen=()
while read -r run rest || [[ -n "${run:-}" ]]; do
  [[ -z "$run" || "$run" == \#* ]] && continue
  if [[ ! "$run" =~ ^[0-9]{1,4}$ || ( -n "$rest" && "$rest" != \#* ) ]]; then
    echo "无效轮次行：$run $rest" >&2
    exit 2
  fi
  printf -v run "%04d" "$((10#$run))"
  if [[ -n "${seen[$run]:-}" ]]; then
    echo "重复轮次：$run" >&2
    exit 2
  fi
  seen[$run]=1
  for candidate in narrow wide; do
    file="$output_dir/data${run}_build_${candidate}.root"
    if [[ -e "$file" || -e "$file.partial" ]]; then
      echo "已有输出，未启动批次：$file；试跑请指定独立输出目录。" >&2
      exit 2
    fi
  done
  if [[ -e "$output_dir/data${run}_build_cpp_complete.json" ]]; then
    echo "已有完成回执，未启动批次：$run" >&2
    exit 2
  fi
  runs+=("$run")
done < "$run_list"
[[ ${#runs[@]} -gt 0 ]] || { echo "轮次名单为空" >&2; exit 2; }

log_dir="$script_dir/terminal_cout/build_$(date +%Y%m%d_%H%M%S)_$$"
mkdir -p -- "$log_dir" "$output_dir" || exit 2
declare -A active=()
failures=0
collect_finished() {
  local pid code
  for pid in "${!active[@]}"; do
    if ! kill -0 "$pid" 2>/dev/null; then
      if wait "$pid"; then code=0; else code=$?; fi
      if (( code != 0 )); then failures=$((failures + 1)); fi
      echo "run ${active[$pid]} finished: exit=$code; log=$log_dir/${active[$pid]}.log"
      unset 'active[$pid]'
    fi
  done
}
for run in "${runs[@]}"; do
  while (( ${#active[@]} >= max_jobs )); do
    collect_finished
    (( ${#active[@]} >= max_jobs )) && sleep 0.2
  done
  "$script_dir/build" "$run" --project "$project_dir" --outdir "$output_dir" >"$log_dir/$run.log" 2>&1 &
  active[$!]="$run"
done
while (( ${#active[@]} > 0 )); do
  collect_finished
  (( ${#active[@]} > 0 )) && sleep 0.2
done
echo "批次结束：${#runs[@]}轮，失败${failures}轮。日志：$log_dir"
(( failures == 0 ))

