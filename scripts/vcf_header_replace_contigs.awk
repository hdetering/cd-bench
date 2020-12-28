#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Replace ##contig entries in VCF file header.
# Purpose: Harmonize sorting order with reference 
# (e.g. when using GATK SortVcf).
#
# usage: awk -f vcf_header_replace_contigs.awk contigs.txt variants.vcf
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-05-21
#------------------------------------------------------------------------------

# read new contig headers
NR==FNR {
  a[++n] = $0;
  next;
}
# check headers in VCF
/^##contig=/ {
  if (c) { # we have encountered a ##contig header before
    next;
  } 
  else { # this is the first ##contig header
    c = 1;
    for (i=1; i<=length(a); i++) {
      print a[i];
    }
    next;
  }
}
# print all other lines
1
