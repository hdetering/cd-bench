#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert CSV with read counts to PyClone-VI input TSV file.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-20
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import vcf

# output header
hdr = "mutation_id sample_id ref_counts alt_counts normal_cn major_cn minor_cn".split()

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone-VI input file from CSV.')
  parser.add_argument('--csv', required=True, type=argparse.FileType(), help='CSV file with somatic mutations and read counts.')
  parser.add_argument('--normal', help='Normal sample id. (ignored for output)')

  args = parser.parse_args()
  return args

def main(args):
  # print output header
  print('\t'.join(hdr))

  # skip input header
  hdr_in = args.csv.readline()
  # parse variant data
  for line in args.csv:
    cols = line.strip().split(',')
    chrom, pos, smp, ref, alt, rc_tot, rc_alt = cols[:7]
    cn_nrm = "2"
    cn_maj = "1"
    cn_min = "1"
    if smp != args.normal:
      id_var = '{}_{}'.format(chrom, pos)
      rc_ref = str(int(rc_tot) - int(rc_alt))
      out = [id_var, smp, rc_ref, rc_alt, cn_nrm, cn_maj, cn_min]
      print('\t'.join(out))

if __name__ == '__main__':
  args = parse_args()
  main(args)
