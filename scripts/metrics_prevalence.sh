#!/bin/bash
# vim: syntax=sh tabstop=2 expandtab
# coding: utf-8

DATADIR="$1"
TOOLDIR="$2"

for rep in $(find "${DATADIR}" -mindepth 1 -maxdepth 1 -type d | sort); do 
  echo $rep
  python3 scripts/metrics_prevalence.py \
    ${rep}/sim/true.clusters.csv ${rep}/sim/true.snvs.csv \
    ${rep}/${TOOLDIR}/inf.clusters.csv ${rep}/${TOOLDIR}/inf.snvs.csv \
  | tee ${rep}/${TOOLDIR}/metrics_prevalence.yml
done

<<'COMMENT'
for rep in $(lsd data/100 | sort | head -1); do 
  python3 scripts/metrics_prevalence.py \
    $rep/sim/true.clusters.csv $rep/sim/true.snvs.csv \
    $rep/phylowgs_cn2_bcftools/inf.clusters.csv $rep/phylowgs_cn2_bcftools/inf.snvs.csv \
  | tee $rep/phylowgs_cn2_bcftools/metrics_prevalence.yml
done

for rep in $(lsd data/100 | sort | head -1); do 
  python3 scripts/metrics_prevalence.py \
    $rep/sim/true.clusters.csv $rep/sim/true.snvs.csv \
    $rep/cloe/cloe.clusters.csv $rep/cloe/cloe.snvs.csv \
  | tee $rep/cloe/metrics_prevalence.yml
done
COMMENT
