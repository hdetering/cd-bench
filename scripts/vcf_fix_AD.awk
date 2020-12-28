#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Replace '.' characters in AD values with '0'.
# Also, add AD value for REF allele if not present.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-07-26
#------------------------------------------------------------------------------

# join an array into a string
function join(array, start, end, sep, result, i)
{
  if (sep == "")
    sep = " "
  else if (sep == SUBSEP) # magic value
    sep = ""
  result = array[start]
  for (i = start + 1; i <= end; i++)
    result = result sep array[i]
  return result
}

BEGIN {
  OFS = "\t";
}

# print header lines
/^#/ {
  print; next;
}

# non-header lines
{
  # split ALT alleles into array
  split($5, alt, ",");
  # split FORMAT tags into array
  split($9, t, ":");
  # loop over samples
  for (i=10; i<=NF; i++) {
    rc_tot = 0; rc_ref = 0; rc_alt = 0;
    # split FORMAT values into array
    split($i, v, ":");
    # loop over FORMAT values
    for (j=1; j<=length(t); j++) {
      if (t[j]=="DP")
        rc_tot = v[j];
      if (t[j]=="AD") {
        if (v[j] != ".")
          gsub(/\./, "0", v[j]);
        # check if all alleles present in AD field
        split(v[j], ad, ",");
        if (length(ad) == length(alt) && rc_tot != ".") { # REF allele missing in AD
          for (k=1; k<=length(ad); k++)
            rc_alt += ad[k];
          rc_ref = rc_tot - rc_alt;
          v[j] = rc_ref","v[j];
        }
      }
    }
    # put FORMAT values back together
    $i = join(v, 1, length(v), ":");
  }
  print;
}
