#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# For each sample in a multi-sample VCF print a '1' if GT==0/1, '0' otherwise.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-13
#------------------------------------------------------------------------------

!/^#/ {
  for (i=10; i<=NF; i++) {
    if ($i ~ /^0\/1/)
      printf 1;
    else 
      printf 0;
  }
  printf "\n";
}
