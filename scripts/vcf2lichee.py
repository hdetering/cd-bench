#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert infos from multiple VCF files to LICHeE input file.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-04-21
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
#import csv
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone YAML input file from VCF and copy-number BED.')
  parser.add_argument('--vcf', required=True, type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP"')
  parser.add_argument('--out', required=True, type=argparse.FileType('wt'), help='Output TSV file.')
  parser.add_argument('--bed', type=argparse.FileType(), help='BED file with allele-specific copy number.')
  parser.add_argument('--use-info', action='store_true', help='Extract read counts from INFO field (overrides --sample)')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  return args

def parse_vcf_record(rec):
  '''
  Parse a multi-sample VCF record.
  Determine major ALT allele and calculate VAF for each sample.
  
  OUTPUT: 
    list with the following columns:
      #chr position description Normal
    empty list if no consensus var.
  '''
  out = []
  # determine major ALT allele
  #import pdb; pdb.set_trace()
  base_cnt = {str(b): 0 for b in rec.ALT}
  for call in rec.samples:
    rc = getattr(call.data, 'AD')
    if rc is None:
      continue
    for i in range(len(rec.ALT)):
      base_cnt[str(rec.ALT[i])] += int(rc[i])
  major_allele = max(base_cnt, key=base_cnt.get)

  # get major allele count for each sample
  idx_maj = ref.ALT.index(major_allele)
  rc = [getattr(c, 'AD')[idx_maj] for c in rec.samples if getattr(c, 'AD') else None]
  # get depth for each sample
  dp = [getattr(c, 'DP') for c in rec.samples]
  # calculate VAFs
  vaf = [float(a)/float(d) for a, d in zip(rc, dp) if (a and d) else 0.0]

  desc = '{}/{}'.format(rec.REF, major_allele)
  out = [ rec.CHROM, str(rec.POS), desc, '0.0' ] + [str(x) for x in vaf]
  return out

def main(args):

  # get samples from input VCF
  rdr = vcf.Reader(filename=args.vcf.name)
  samples = rdr.samples
  
  # write header
  fh_out = args.out
  hdr = "#chr position description Normal".split() + samples
  fh_out.write('\t'.join(hdr) + '\n')

  # parse variants
  for rec in rdr:
    var_data = parse_vcf_record(rec)
    if True: #len(var_line) == len(hdr):
      sys.stdout.write('\t'.join(var_data) + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args)
