#!/bin/bash
# vim: syntax=sh tabstop=2 shiftwidth=2 expandtab
# coding: utf-8

DO_CLUST=true
DO_PREV=true
DO_PHYLO=true
SCRIPTS="/mnt/netapp2/posadalab2/uvibehde/cd-bench/scripts"
#DATADIR="data/100"
DATADIR="$1"
#ANADIR="analysis/100"
ANADIR="$2"

tools_clust=(
  "cloe_cn2_bcftools"
  "cloe_cn2_seqz"
  "clonefinder_cn2_bcftools"
  "clonefinder_cn2_seqz"
  "lichee_cn2_bcftools"
  "lichee_cn2_seqz"
  "phylowgs_cn2_bcftools"
  "phylowgs_cn2_seqz"
  "pyclone_cn2_bcftools"
  "pyclone_cn2_seqz"
  "pyclonevi_cn2_bcftools"
  "pyclonevi_cn2_seqz"
  "pyclonevi_cn"
  "sciclone_cn2_bcftools"
  "sciclone_cn2_seqz"
)
tools_prev=("${tools_clust[@]}")
tools_phylo=(
  "cloe_cn2_bcftools"
  "cloe_cn2_seqz"
  "clonefinder_cn2_bcftools"
  "clonefinder_cn2_seqz"
  "lichee_cn2_bcftools"
  "lichee_cn2_seqz"
  "phylowgs_cn2_bcftools"
)
tools_phylo=(
  "cloe_cn2_bcftools"
  "cloe_cn2_seqz"
  "clonefinder_cn2_bcftools"
  "clonefinder_cn2_seqz"
  "lichee_cn2_bcftools"
  "lichee_cn2_seqz"
)

if [ "${DO_CLUST}" = true ]; then
  for tool in ${tools_clust[@]}; do
    bash ${SCRIPTS}/gather_clustering.sh ${DATADIR} ${tool} ${ANADIR}
  done
fi

if [ "${DO_PREV}" = true ]; then
  for tool in ${tools_prev[@]}; do
    bash ${SCRIPTS}/gather_prevalence.sh ${DATADIR} ${tool} ${ANADIR}
  done
fi
if [ "${DO_PHYLO}" = true ]; then
  for tool in ${tools_phylo[@]}; do
    bash ${SCRIPTS}/gather_phylo.sh ${DATADIR} ${tool} ${ANADIR}
  done
fi
