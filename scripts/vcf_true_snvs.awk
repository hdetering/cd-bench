#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Extract true SNVs and their clone assignments from VCF and mutation matrix.
#  usage: awk -f vcf_true_snvs.awk somatic.vcf mm.csv
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-21
#------------------------------------------------------------------------------

BEGIN {
  OFS = ",";
}
# parse true SNVs
NR == FNR {
  if ($0 !~ /^#/) {
    
  }
  next;
}
# parse mutation matrix
