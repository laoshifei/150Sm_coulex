#!/usr/bin/bash
# Legacy prefix-sample workflow; prefer calibrate and explicit batch_new.sh inputs.

max_jobs=${MAX_JOBS:-1}

while read run; do
  mkdir -p fig/$(printf "%04d" $run) 
  root -b -q "tree2hist_check_fast.cpp($run)" &

  while (( $(jobs -p | wc -l) >= $max_jobs )); do
    wait -n
  done
done < "../pars/run.txt"

