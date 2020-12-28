#!/bin/bash
# vim: syntax=sh tabstop=2 shiftwidth=2 expandtab
# coding: utf-8

# INPUT:
# 1st param: clusters file
# sample_id	cluster_id size mean std
IN_CLUST=$1
# 2nd param: loci file
# mutation_id sample_id cluster_id cellular_prevalence cellular_prevalence_std variant_allele_frequency
IN_LOCI=$2

# OUTPUT:
# output goes to results directory
OUTDIR=$(dirname ${IN_CLUST})
# 1st output is clusters file
# id_cluster,id_sample,freq
OUT_CLUST="${OUTDIR}/inf.clusters.csv"
# 2nd output is snvs file
# id_mut,id_cluster
OUT_SNV="${OUTDIR}/inf.snvs.csv"

echo "-- INPUT --"
echo ${IN_CLUST}
echo ${IN_LOCI}
echo "-- OUTPUT --"
echo ${OUT_CLUST}
echo ${OUT_SNV}

# convert clusters file
tail -n+2 ${IN_CLUST} |
awk '
BEGIN { 
  OFS = ",";
  print "id_cluster,id_sample,freq";
}
{ print $2,$1,$4; }
' > ${OUT_CLUST}

# convert loci file
tail -n+2 ${IN_LOCI} |
awk '
BEGIN { 
  OFS = ",";
  print "chrom_pos,id_cluster";
}
{ print $1,$3 }
' > ${OUT_SNV}
