#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Calculate BAF and LRR from a paired normal-tumor VCF.
#
# Trim ALT alleles to keep allele with max AD (random choice to break ties). 
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-05-17
#------------------------------------------------------------------------------

import os, sys
import argparse
import math
import random
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Calculate BAF and LRR from a paired normal-tumor VCF.')
  parser.add_argument('vcf', type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP"')
  parser.add_argument('normal', help='Normal sample id.')

  args = parser.parse_args()
  return args

def calc_baf(rec):
  '''
  Calculate B-allele frequency (BAF) from multi-sample VCF record.
  
  In case of multiple ALT alleles, pick allele with the highest AD overall.
  Break ties in AD by randomly choosing an allele.

  Returns:
    - list with BAF for each sample
  '''

  # sum up AD values for each ALT allele
  sum_ad = [0] * len(rec.ALT)
  for call in rec.samples:
    for i in range(len(rec.ALT)):
      if call.data.AD and call.data.AD[i+1]:
        sum_ad[i] += call.data.AD[i+1]
  # determine major ALT allele
  max_ad = max(sum_ad)
  i_max = [i for i,v in enumerate(sum_ad) if v == max_ad]
  i_maj = i_max[0] if len(i_max)==1 else random.choice(i_max)
  # calculate BAF of major allele for all samples
  #import pdb; pdb.set_trace()
  baf = [float(c.data.AD[i_maj+1])/c.data.DP if c.data.AD else None for c in rec.samples]
  
  return baf
  
def parse_record(rec, idx_normal, mean_dp):
  '''
  Parse multi-sample VCF record.
  
  Apply the following filters:
    - Normal sample DP must not be Null
  Do these modifications:
    - keep only ALT allele with highest AD
    - set BAF to AD(ALT)/DP
    - set LRR to log2(DP/mean_dp)
  '''
  # only proceed if normal sample has variant
  if rec.samples[idx_normal].data.DP is None:
    return None
  # create new CallData class (NamedTuple)
  old_keys = list(rec.samples[idx_normal].data._fields)
  cd = vcf.parser.make_calldata_tuple(old_keys + ['BAF', 'LRR'])
  # update VCF record
  rec.FORMAT += ':BAF:LRR'
  baf = calc_baf(rec)
  for i in range(len(rec.samples)):
    # calculate LRR from DP
    lrr = None
    if not rec.samples[i].data.DP is None:
      rr = float(rec.samples[i].data.DP) / mean_dp
      lrr = math.log2(rr) if rr > 0 else 0.0
    s_lrr = '{:0.4f}'.format(lrr) if lrr else None
    s_baf = '{:0.4f}'.format(baf[i]) if baf[i] else None
    old_vals = list(rec.samples[i].data._asdict().values())
    data = cd(*(old_vals + [s_baf, s_lrr]))
    rec.samples[i].data = data

  return rec

def main(args):
  '''
  Parse multi-sample VCF file, per chromosome do 2 passes:
    1. calculate mean read depth of normal sample
    2. calculate BAF and LRR for major ALT allele
  '''
  rdr = vcf.Reader(filename=args.vcf.name)
  idx_normal = rdr.samples.index(args.normal) # ValueError if not present

  # extend VCF header with additional FORMAT fields
  hdr = rdr
  hdr.formats['BAF'] = vcf.parser._Format('BAF', 1, 'Float', 'B-allele frequency')
  hdr.formats['LRR'] = vcf.parser._Format('LRR', 1, 'Float', 'Log response ratio (log2 of ratio (DP/(chrom-wise mean DP_normal))') 

  wtr = vcf.Writer(sys.stdout, hdr)
  # 1st pass
  recs = []
  chrom = ""
  sum_dp = 0
  num_rec = 0
  mean_dp = 0.0
  for rec in rdr:
    if rec.CHROM == chrom:
      if rec.samples[idx_normal].data.DP:
        recs.append(rec)
        sum_dp += rec.samples[idx_normal].data.DP
        num_rec += 1
    else: # passed end of chromosome
      if len(recs) > 0:
        mean_dp = float(sum_dp) / num_rec
        # 2nd pass
        for r in recs:
          rec_out = parse_record(r, idx_normal, mean_dp)
          if rec_out:
            wtr.write_record(rec_out)
      # prepare for processing next chrom's variants
      recs = [rec]
      chrom = rec.CHROM
      if rec.samples[idx_normal].data.DP:
        sum_dp = rec.samples[idx_normal].data.DP
        num_rec = 1
      else:
        sum_dp = 0
        num_rec = 0

if __name__ == '__main__':
  args = parse_args()
  main(args)
