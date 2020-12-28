#!/usr/bin/env python2
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert multi-sample VCF to SPRUCE input file.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-16
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import vcf
import spruce_convert

out_hdr = '#id_sample lbl_sample id_mut lbl_mut vaf_lb vaf_mean vaf_ub x y mu'

def parse_args():
  parser = argparse.ArgumentParser(description='Create SPRUCE input file from VCF.')
  parser.add_argument('--vcf', required=True, type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP"')
  parser.add_argument('--out', required=True, type=argparse.FileType('wt'), help='Output TSV file.')
  parser.add_argument('--normal', help='Normal sample id. (Will be ignored if present)')

  args = parser.parse_args()
  return args

def parse_vcf_record(rec, lbl_samples):
  '''
  Parse a multi-sample VCF record.
  Determine major ALT allele and calculate VAF for each sample.
  
  OUTPUT: 
    list with the following elements:
      CHROM_POS REF ALT (S1:ref, S1:alt) [(S2:ref, S2:alt) [...]]
    empty list if no consensus var.
  '''
  out = []
  # reorder samples
  idx_smp = [rec._sample_indexes[x] for x in lbl_samples]
  samples = [rec.samples[i] for i in idx_smp if i < len(rec.samples)]
  # determine major ALT allele
  base_cnt = {str(b): 0 for b in rec.ALT}
  #for call in samples:
  for i in idx_smp:
    if i >= len(rec.samples):
      continue
    rc = getattr(rec.samples[i].data, 'AD') if hasattr(rec.samples[i].data, 'AD') else None
    if rc is None:
      continue
    for i in range(len(rec.ALT)):
      base_cnt[str(rec.ALT[i])] += int(rc[i+1] or 0)
  major_allele = max(base_cnt, key=base_cnt.get)

  # get major allele count for each sample
  idx_maj = rec.ALT.index(major_allele)
  rc = [None] * len(idx_smp)
  dp = [None] * len(idx_smp)
  for i in idx_smp:
    if i < len(rec.samples):
      cdata = rec.samples[i].data
      ad = getattr(cdata, 'AD') if hasattr(cdata, 'AD') else None
      if ad:
        rc[i] = ad[idx_maj+1]
      dp[i] = getattr(cdata, 'DP')

  # compile REF/ALT counts
  alt = [  x if x else 0 for   x in rc]
  ref = [d-a if d else 0 for a,d in zip(alt, dp)]

  id_var = '{}_{}'.format(rec.CHROM, str(rec.POS))
  out = [id_var] + [(r, a) for r,a in zip(ref, alt)]
  return out

def main(args):
  # get samples from input VCF
  rdr = vcf.Reader(filename=args.vcf.name)
  samples = rdr.samples
  assert args.normal in samples
  # make sure normal sample is first in list
  if args.normal:
    samples = [args.normal] + [x for x in samples if x != args.normal]

  # parse variants
  mut_rc = [] # [id_mut, (ref_RN, alt_RN), (ref_R1, alt_R1), ...]
  for rec in rdr:
    rc_data = parse_vcf_record(rec, samples)
    # make sure that mutation is present in at least one sample
    rc_alt = sum([a for r,a in rc_data[2:]])
    if rc_alt > 0:
      mut_rc.append(rc_data)

  # construct header
  fh_out = args.out
  hdr  = '{} # number of samples\n'.format(len(samples)-1)
  hdr += '{} # number of SNVs\n'.format(len(mut_rc))
  hdr += '\t'.join(out_hdr.split()) + '\n'
  fh_out.write(hdr)
  # write variants
  for i, mut in enumerate(mut_rc):
    for j in range(len(samples)-1):
      rc_ref, rc_alt = mut[j+2]
      # calculate limits of confidence interval (highest posterior density region)
      mode = 0; lb = 0; ub = 0
      if rc_alt+rc_ref > 0:
        mode, lb, ub = spruce_convert.binomial_hpdr(rc_alt, rc_alt+rc_ref, 0.999)
      out  = '{}\t{}'.format(j, samples[j+1])
      out += '\t{}\t{}'.format(i, mut[0])
      out += '\t{:.4f}\t{:.4f}\t{:.4f}'.format(lb, mode, ub)
      # CN1, CN2, prev
      out += '\t{}\t{}\t{}'.format(1, 1, 1)
      fh_out.write(out + '\n')

if __name__ == '__main__':
  args = parse_args()

  main(args)
