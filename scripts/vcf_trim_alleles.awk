#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# For loci with multiple ALT alleles, select the one with highest AD.
# Ties are resolved by random selection.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-05-14
#------------------------------------------------------------------------------

BEGIN {
  srand(systime()); # set random seed
  OFS = "\t";
}
# non-header lines
!/^#/ {
  # get ALT alleles
  split($5, alt, ",");
  if (length(alt)>1) {
    # get FORMAT fields
    split($9, tags, ":");
    # get sample values
    split($10, vals, ":");
    # parse FORMAT tags
    for (i in tags) {
      if (tags[i] == "AD") {
        split(vals[i], ad, ",");
        ad_new = ad[1]; # first AD value is REF allele
        # determine max AD value
        m = 0;
        for (j=2; j<=length(ad); j++) {
          if (ad[j] > m) m = ad[j];
        }
        # determine allele(s) having max AD
        n = 0; split("", a); # clear array "a"
        for (j=1; j<length(ad); j++) {
          if (ad[j+1] == m) a[++n] = j;
        }
        # pick random allele index from list (to break ties)
        r = int(rand()*length(a)) + 1;
        idx_m = a[r];
        # new values to output
        vals[i] = ad[1]","m;
        $5 = alt[idx_m];
      }
      if (i==1) new_vals = vals[i];
      else new_vals = new_vals":"vals[i];
    }
    $10 = new_vals;
  }
}
# print line (modified or not)
{ print; }
