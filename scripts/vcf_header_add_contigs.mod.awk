#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Add ##contig entries in VCF file header.
# Headers will be inserted before #CHROM header line-
# Purpose: bcftools annotate requires that ##contig headers be present
#
# usage: awk -f vcf_header_add_contigs.awk contigs.txt variants.vcf
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
# omit faulty AD format tag specification (include correct one in contigs.txt)
/^##FORMAT=<ID=AD,Number=A/ {
  next;
}
/^#CHROM/ {
  for (i=1; i<=length(a); i++) {
    print a[i];
  }
}
# print all other lines
1
