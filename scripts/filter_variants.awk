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
# modified : 2020-03-18
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
  split($8, tags, ";");
  for (i in tags) {
    split(tags[i], kv, "=");
    if (kv[1] == "DP") {
      rc_tot = kv[2];
    }
    else if (kv[1] == "AC") {
      rc_alt = kv[2];
    }
  }

  # print variant line if filter criteria met
  if (rc_alt>=rc_min && rc_alt/rc_tot>=vaf_min) {
    print;
  }
}
