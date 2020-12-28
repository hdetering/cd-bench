#!/bin/bash
DATADIR=$1

for rep in $(find ${DATADIR} -mindepth 1 -maxdepth 1 -type d | sort); do
  n_mut=$(cat ${rep}/log/lichee_cn2_strict.log | grep 'somatic SNVs' | cut -d" " -f7)
  n_tree=$(grep 'valid tree' ${rep}/log/lichee_cn2_strict.log | tail -1 | awk '{print $2}')
  echo "$(basename ${rep}),${n_mut},${n_tree}"; done
