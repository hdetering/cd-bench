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
FILENAME=$2

for rep in $(find ${DATADIR} -mindepth 1 -maxdepth 1 -type d | sort); do
  id_rep=$(basename ${rep})
  fn="${rep}/sim/${FILENAME}"
  # check if header was output before
  if [ -z "${hdr+x}" ]; then
    hdr=$(head -1 ${fn})
    echo "id_rep,${hdr}"
  fi
  if [[ -f "${fn}" ]]; then
    tail -n+2 ${fn} | awk -v FS="," -v rep=${id_rep} '{print rep","$0}'
  fi
done
