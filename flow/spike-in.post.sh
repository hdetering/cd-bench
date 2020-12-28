#!/bin/bash
# vim: syntax=sh tabstop=2 shiftwidth=2 expandtab
# coding: utf-8

# run these steps after executing
#   simtools/simtools.py spikein ...

SCRIPTS="/mnt/netapp2/posadalab2/uvibehde/cd-bench/scripts"
wd="/mnt/netapp2/posadalab2/uvibehde/cd-bench/data/200/33adb830-fbf3-11ea-abcf-a0000220fe80"
wd=$1
cd ${wd}/sim/bam

module load bcftools
set -x

# make sure that VCFs are
# - bgzip'ed (gzip'ed won't work)
# - bcftools sort'ed
# - bcftools index'ed
for f in R*.rc.vcf.gz; do
  if [ $(htsfile $f | grep -c "BGZF") -eq 0 ]; then
    gunzip $f
    bgzip ${f%.gz}
  fi
  bcftools index -f $f
done
for f in R*.rc.spikein.vcf.gz; do
  if [ $(htsfile $f | grep -c "BGZF") -eq 0 ]; then
    mv $f ${f%.gz}
    bgzip ${f%.gz}
  fi
  bcftools sort -Oz $f > ${f%.vcf.gz}.sort.vcf.gz
  bcftools index ${f%.vcf.gz}.sort.vcf.gz
done

# merge VCFs
for f in R*.rc.vcf.gz; do
  bcftools merge --force-samples $f ${f%.vcf.gz}.spikein.sort.vcf.gz \
  | awk -f ${SCRIPTS}/vcf_merge_single.awk \
  | bgzip -c \
  > ${f%.rc.vcf.gz}_tail.rc.vcf.gz
done

# filter variants
for f in R*_tail.rc.vcf.gz; do
  zcat $f | awk -f ${SCRIPTS}/filter_variants.awk \
  | bgzip -c \
  > ${f%.vcf.gz}.filt.vcf.gz
  bcftools index ${f%.vcf.gz}.filt.vcf.gz 
done

# create links for CN beds
cd -
cd ${wd}/sim/bed
for f in R*.cn.bed.gz; do 
  ln -s $f ${f%.cn.bed.gz}_tail.cn.bed.gz
  ln -s ${f}.tbi ${f%.cn.bed.gz}_tail.cn.bed.gz.tbi
done
