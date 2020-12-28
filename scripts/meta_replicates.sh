#!/usr/bin/bash
# vim: syntax=sh tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Collect meta data about replicates.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-12-13
#------------------------------------------------------------------------------

# directory containing replicate dirs can be passed, otherwise use cwd.
#DATADIR=${1:-"."}

echo "#id_rep,n_clones,n_samples,purity,depth,n_muts,cnv_rate,n_snvs"
for cohort in data/???; do
for rep in $(find ${cohort} -mindepth 1 -maxdepth 1 -type d | sort); do
  id_rep=$(basename ${rep})
  n_clones=$(grep "nclones" ${rep}/meta.yml | awk '{print $2}')
  n_samples=$(grep "nsamples" ${rep}/meta.yml | awk '{print $2}')
  purity=$(grep "purity" ${rep}/meta.yml | awk '{print $2}')
  depth=$(grep "seq-coverage" ${rep}/config.yml | awk '{print $2}')
  muts=$(grep "mut-som-num" ${rep}/config.yml | awk '{print $2}')
  cnvs=$(grep "mut-som-cnv-ratio" ${rep}/config.yml | awk '{print $2}')
  snvs=$(grep -cv "^#" ${rep}/sim/somatic.vcf)
  echo "${id_rep},${n_clones},${n_samples},${purity},${depth},${muts},${cnvs},${snvs}"
done
done
