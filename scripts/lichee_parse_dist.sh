#!/bin/bash
# vim: syntax=sh tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
DATADIR=$1
SUBDIR="lichee/cn2_relaxed"

# write header
echo "id_rep,caset_union,caset_isect,disc_union,disc_isect"
for rep in $(find ${DATADIR} -mindepth 1 -maxdepth 1 -type d | sort); do
  caset_union=""
  caset_isect=""
  disc_union=""
  disc_isect=""
  if [[ -f "$rep/${SUBDIR}/snv.vaf.tsv.trees.txt" ]]; then
    caset_union=$(tail -1 $rep/${SUBDIR}/lichee.tree.0.CASet_union.dist)
    caset_isect=$(tail -1 $rep/${SUBDIR}/lichee.tree.0.CASet_isect.dist)
    disc_union=$(tail -1 $rep/${SUBDIR}/lichee.tree.0.DISC_union.dist)
    disc_isect=$(tail -1 $rep/${SUBDIR}/lichee.tree.0.DISC_isect.dist)
  fi
  echo "$(basename ${rep}),${caset_union},${caset_isect},${disc_union},${disc_isect}"; done
