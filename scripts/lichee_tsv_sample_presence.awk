#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# For each sample in a LICHeE input TSV print a '1' if VAF>0, '0' otherwise.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-16
#------------------------------------------------------------------------------

BEGIN {
  OFS = "\t";
}

!/^#/ {
  pres = "";
  for (i=4; i<=NF; i++) {
    if ($i > 0.0)
      pres = pres 1;
    else 
      pres = pres 0;
  }
  print pres,$0;
}
