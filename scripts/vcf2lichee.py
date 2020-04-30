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
  #parser.add_argument('--bed', type=argparse.FileType(), help='BED file with allele-specific copy number.')
  parser.add_argument('--normal', help='Normal sample id. (default: assume VAF=0.0 for Normal)')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  return args

def parse_vcf_record(rec, lbl_samples, add_normal):
  '''
  Parse a multi-sample VCF record.
  Determine major ALT allele and calculate VAF for each sample.
  
  OUTPUT: 
    list with the following columns:
      #chr position description Normal
    empty list if no consensus var.
  '''
  out = []
  # reorder samples
  samples = [rec.samples[rec._sample_indexes[x]] for x in lbl_samples]
  # determine major ALT allele
  base_cnt = {str(b): 0 for b in rec.ALT}
  for call in samples:
    rc = getattr(call.data, 'AD')
    if rc is None:
      continue
    for i in range(len(rec.ALT)):
      base_cnt[str(rec.ALT[i])] += int(rc[i])
  major_allele = max(base_cnt, key=base_cnt.get)

  # get major allele count for each sample
  idx_maj = rec.ALT.index(major_allele)
  #import pdb; pdb.set_trace()
  rc = [getattr(c.data, 'AD')[idx_maj] if getattr(c.data, 'AD') else None for c in samples]
  # get depth for each sample
  dp = [getattr(c.data, 'DP') for c in samples]
  # calculate VAFs
  vaf = [float(a)/float(d) if (a and d) else 0.0 for a, d in zip(rc, dp)]
  if add_normal:
    vaf = [0.0] + vaf

  desc = rec.ID if rec.ID else '.' 
  #desc = '{}/{}'.format(rec.REF, major_allele)
  out = [rec.CHROM, str(rec.POS), desc] + ['{:.4f}'.format(x) for x in vaf]
  return out

def main(args):

  # get samples from input VCF
  rdr = vcf.Reader(filename=args.vcf.name)
  samples = rdr.samples
  do_add_normal = (args.normal is None)
  assert do_add_normal or args.normal in samples
  # make sure normal sample is first in list
  if args.normal:
    samples = [args.normal] + [x for x in samples if x != args.normal]

  # write header
  fh_out = args.out
  hdr = "#chr position description Normal".split() + samples
  fh_out.write('\t'.join(hdr) + '\n')

  # parse variants
  for rec in rdr:
    var_data = parse_vcf_record(rec, samples, do_add_normal)
    if True: #len(var_line) == len(hdr):
      fh_out.write('\t'.join(var_data) + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args)
