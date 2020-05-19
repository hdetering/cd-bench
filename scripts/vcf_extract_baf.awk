#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Extract B-allele frequency (BAF) for biallelic variants
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-05-14
#------------------------------------------------------------------------------

# non-header lines
!/^#/ {
  # get ALT alleles
  split($5, alt, ",");
  if (length(alt) == 1) { # can determine BAF for biallelic SNVs only
    # get FORMAT fields
    split($9, tags, ":");
    # get sample values
    split($10, vals, ":");
    # parse FORMAT tags
    for (i in tags) {
      if (tags[i] == "AD") {
        split(vals[i], ad, ",");
        baf = ad[2]/(ad[1]+ad[2]);
        printf("%s\t%d\t%.04f\n", $1, $2, baf);
        next;
      }
    }
  }
}
