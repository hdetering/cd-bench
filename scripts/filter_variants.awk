#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Filter variants in a given VCF file by the following criteria:
#  - ALT allele read count >= minimum
#  - VAF >= minimum
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-05-18
#------------------------------------------------------------------------------

BEGIN {
  rc_min = 3
  vaf_min = 0.05
}
# header lines
/^#/ {
  print;
  next;
}
# non-header lines
{
  rc_tot = 0;
  rc_alt = 0;
  # split info fields
  split( $9, tags, ":");
  split($10, vals, ":");
  for (i in tags) {
    if (tags[i] == "DP") {
      rc_tot = vals[i];
    }
    else if (tags[i] == "AD") {
      split(vals[i], ad, ",");
      # get max AD among ALT alleles
      max_ad = 0;
      for (i=2; i<=length(ad); i++) {
        if (ad[i] > max_ad) max_ad = ad[i];
      }
      rc_alt = max_ad;
    }
  }
  # print variant line if filter criteria met
  if (rc_alt>=rc_min && rc_alt/rc_tot>=vaf_min) {
    print;
  }
}
