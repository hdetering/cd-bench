#!/bin/bash
# vim: syntax=sh tabstop=2 expandtab
# coding: utf-8

DATADIR="$1"
TOOLDIR="$2"

for rep in $(find "${DATADIR}" -mindepth 1 -maxdepth 1 -type d | sort); do 
  python3 scripts/metrics_phylogeny.py \
    --CASet_isect --CASet_union --DISC_isect --DISC_union \
    --collapse --ignore-homoplasy \
    ${rep}/sim/true.tree.csv ${rep}/sim/true.snvs.csv \
    ${rep}/${TOOLDIR}/inf.trees.csv ${rep}/${TOOLDIR}/inf.snvs.csv \
  | tee ${rep}/${TOOLDIR}/metrics_phylogeny.yml
done

<<'COMMENT'
for rep in $(lsd data/100 | sort); do 
  python3 scripts/metrics_phylogeny.py \
    --CASet_isect --CASet_union --DISC_isect --DISC_union \
    --collapse --ignore-homoplasy \
    $rep/sim/true.tree.csv $rep/sim/true.snvs.csv \
    $rep/cloe_cn2_bcftools/inf.trees.csv $rep/cloe_cn2_bcftools/inf.snvs.csv \
  | tee $rep/cloe_cn2_bcftools/metrics_phylogeny.yml
done

for rep in $(lsd data/100 | sort); do 
  python3 scripts/metrics_phylogeny.py \
    --CASet_isect --CASet_union --DISC_isect --DISC_union \
    --collapse --ignore-homoplasy \
    $rep/sim/true.tree.csv $rep/sim/true.snvs.csv \
    $rep/cloe/cloe.trees.csv $rep/cloe/cloe.snvs.csv \
  | tee $rep/cloe/metrics_phylogeny.yml
done
COMMENT
