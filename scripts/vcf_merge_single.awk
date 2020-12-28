#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Merge VCFs to spike in variants.
#
# Collisions are resolved by
#  - summing up AD values for each ALT allele
#  - update DP to sum of ADs
#
# INPUT: 
#  - VCF with two samples columns, actually representing the same sample:
#    1. orginal variants
#    2. new variants to spike in
# OUTPUT:
#  - Merged VCF containing only the first sample, but with updated AD values
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-10-07
#------------------------------------------------------------------------------

# META headers
/^##/ {
  print;
  next;
}
# main header line
/^#CHROM/ {
  for (i=1; i<10; i++) {
    printf("%s\t", $i);
  }
  printf("%s\n", $10);
  next;
}
# variant record
!/^#/ {
  for (i=1; i<10; i++) {
    printf("%s\t", $i);
  }
  if ($11 ~ /^\./) { # spike-in sample doesn't have locus
    # just print original variant
    printf("%s\n", $10);
    next;
  }
  if ($10 ~ /^\./) { # original sample doesn't have locus
    # just print spike-in variant
    printf("%s\n", $11);
    next;
  }
  # if we got here, both original and spike-in sample have a variant
  # need to resolve collision
  split($10, ori, ":");
  split($11, spk, ":");
  # extract AD values for each ALT allele
  split(ori[2], ori_ad, ",");
  split(spk[2], spk_ad, ",");
  # reconstruct FORMAT/AD field, summing across samples
  ad = ori_ad[1];
  for (i=2; i<=length(ori_ad); i++) {
    if (ori_ad[i] == ".") ori_ad[i] = 0;
    if (spk_ad[i] == ".") spk_ad[i] = 0;
    ad = ad","(ori_ad[i] + spk_ad[i]);
  }
  ori[2] = ad;
  # update FORMAT/DP field
  split(ad, ads, ",");
  dp = ads[1];
  for (i=2; i<=length(ads); i++) dp += ads[i];
  fmt = dp;
  for (i=2; i<=length(ori); i++) {
    fmt = fmt":"ori[i];
  }
  printf("%s\n", fmt);
}
