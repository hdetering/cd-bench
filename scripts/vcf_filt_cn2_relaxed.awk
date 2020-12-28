#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Filter multi-sample VCF file for the following criteria:
#  - FORMAT/CNI == 2
# 
# Samples which don't meet the criteria have their FORMAT values reset to '.'
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-02
#------------------------------------------------------------------------------

BEGIN {
  OFS = "\t";
}
# print header lines
/^#/ {
  print;
  next;
}
# non-header lines
{
  # parse FORMAT tags
  split($9, tags, ":");
  for (i=1; i<=length(tags); i++) {
    if (tags[i] == "CNI")
      idx_cni = i;
    else if (tags[i] == "DP")
      idx_dp = i;
  }
  n = 0; # number of samples having variant (before filtering)
  has_var = 0; # flag indicating if current sample has variant
  # can filter be applied?
  if (idx_cni) {
    # check all samples' FORMAT values
    for (i=10; i<=NF; i++) {
      split($i, vals, ":");
      if (vals[idx_dp] > 0)
        has_var = 1;
      # sanity check: values should match tags (otherwise hands off)
      if (length(vals) == length(tags)) {
        if (vals[idx_cni] != 2) { # filter this call out!
          has_var = 0;
          new_vals = ".";
          for (j=2; j<=length(vals); j++)
            new_vals = new_vals ":.";
          $i = new_vals;
        }
      }
    n += has_var;
    }
  }
}
# print any line (modified or not)
{
  if (n>0) print;
}
