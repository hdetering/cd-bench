#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert infos from multisample VCF file to SciClone input files.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-22
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys
import argparse
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Create SciClone CSV input files from VCF.')
  parser.add_argument('--vcf', required=True, type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP"')
  parser.add_argument('--outdir', required=True, help='Output directory.')
  parser.add_argument('--normal', help='Normal sample id. (default: assume VAF=0.0 for Normal)')

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
      
  #rc = [getattr(c.data, 'AD')[idx_maj+1] if getattr(c.data, 'AD') else None for c in samples]
  # get depth for each sample
  #dp = [getattr(c.data, 'DP') for c in samples]
  # calculate VAFs
  vaf = [float(a)/float(d) if (a and d) else 0.0 for a, d in zip(rc, dp)]
  if add_normal:
    vaf = [0.0] + vaf

  desc = rec.ID if rec.ID else '.' 
  #desc = '{}/{}'.format(rec.REF, major_allele)
  out = [rec.CHROM, str(rec.POS), desc] + [(rc[i], dp[i], x) for i, x in enumerate(vaf)]
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

  # parse variants
  num_loc = 0
  sum_depth = 0
  sample_data = {x: [] for x in samples[1:]}
  for rec in rdr:
    var_data = parse_vcf_record(rec, samples, do_add_normal)
    chrom = var_data[0]
    pos = var_data[1]
    for i in range(1, len(samples)):
      alt, tot, vaf = var_data[3+i]
      if not tot is None:
        num_loc += 1
        sum_depth += tot
      var_dat = (chrom, pos, alt, tot, vaf)
      sample_data[samples[i]].append(var_dat)
  mean_depth = round(sum_depth / num_loc)  

  # output header
  hdr = "chr pos refCount varCount VAF".split()
  # write one output file for each sample
  for id_smp in sorted(sample_data.keys()):
    fn_out = os.path.join(args.outdir, '{}.vaf.csv'.format(id_smp))
    with open(fn_out, 'wt') as f:
      f.write(','.join(hdr) + '\n')
      for var_dat in sample_data[id_smp]:
        chrom, pos, alt, tot, vaf = var_dat
        if tot is None: # if variant is not present in this sample, fill in mean depth
          tot = mean_depth
          alt = 0
        out_dat = (chrom, pos, str(tot-alt), str(alt), '{:.4f}'.format(vaf*100))
        f.write(','.join(out_dat) + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args)
