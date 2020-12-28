#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert infos from multi-sample VCF file to CloneFinder input file.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-11-09
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Create CloneFinder input file from VCF.')
  parser.add_argument('--vcf', required=True, type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP"')
  parser.add_argument('--out', required=True, type=argparse.FileType('wt'), help='Output TSV file.')
  #parser.add_argument('--bed', type=argparse.FileType(), help='BED file with allele-specific copy number.')
  parser.add_argument('--normal', help='Normal sample id. (default: assume VAF=0.0 for Normal)')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

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

  desc = rec.ID if rec.ID else '.'
  id_var = '{}_{}'.format(rec.CHROM, str(rec.POS)) 
  out = [id_var, rec.REF, major_allele] + [(r, a) for r,a in zip(ref, alt)]
  return out

def main(args):

  # get samples from input VCF
  rdr = vcf.Reader(filename=args.vcf.name)
  samples = rdr.samples
  # make sure normal sample is first in list
  if args.normal:
    assert args.normal in samples
    samples = [args.normal] + [x for x in samples if x != args.normal]

  # parse variants
  var_data = []
  smp_dp   = [0] * len(samples) # cumulative depth for each sample
  smp_nloc = [0] * len(samples) # number of loci for each sample
  for rec in rdr:
    vdat = parse_vcf_record(rec, samples)
    var_data.append(vdat)
    for i in range(len(samples)):
      dp = sum(vdat[3+i])
      if dp > 0:
        smp_dp[i] += dp
        smp_nloc[i] += 1
  # calculate mean depth for each sample
  smp_mean_dp = [round(dp/n) if n>0 else 0 for dp,n in zip(smp_dp, smp_nloc)]

  # write header
  fh_out = args.out
  hdr = "#SNVID Wild Mut".split()
  for id_sample in samples:
    hdr += ['{}:ref'.format(id_sample), '{}:alt'.format(id_sample)]
  fh_out.write('\t'.join(hdr) + '\n')

  # write variants
  for vdat in var_data:
    line = '\t'.join(vdat[:3])
    for i in range(len(samples)):
      ref, alt = vdat[3+i]
      if ref + alt == 0: # if locus is absent use mean sample depth as REF count
        ref = smp_mean_dp[i]
      line += '\t{}\t{}'.format(ref, alt)
    fh_out.write(line + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args)
