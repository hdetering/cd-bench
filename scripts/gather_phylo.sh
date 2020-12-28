#!/bin/bash
# vim: syntax=sh tabstop=2 shiftwidth=2 expandtab
# coding: utf-8

if [[ $# -lt 2 ]]; then
  echo
  echo "usage: bash $0 DATADIR OUTDIR"
  echo
  exit 1
fi

DATADIR=$1
TOOL=$2
OUTDIR=$3
SCRIPTS="$HOME/code/cd-bench/scripts"
OUT_YML="${OUTDIR}/metrics_phylo.${TOOL}.yml"
OUT_CSV="${OUTDIR}/metrics_phylo.${TOOL}.csv"

for rep in $(find ${DATADIR} -mindepth 1 -maxdepth 1 -type d | sort); do
  id_rep=$(basename ${rep})
  fn="${rep}/${TOOL}/metrics_phylogeny.yml"
  if [[ -f "${fn}" ]]; then
    echo "- {"
    echo "id_rep: ${id_rep},"
    echo "tool: ${TOOL},"
    cat ${fn} | awk '{print $0","}'
    echo "}"
  fi
done > ${OUT_YML}

python3 ${SCRIPTS}/yaml2csv.py ${OUT_YML} > ${OUT_CSV}
