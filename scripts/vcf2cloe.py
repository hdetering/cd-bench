#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert infos from multisample VCF file to Cloe input files.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-30
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys
import argparse
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Create Cloe input files from VCF.')
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
      
  # calculate VAFs
  vaf = [float(a)/float(d) if (a and d) else 0.0 for a, d in zip(rc, dp)]

  out = [rec.CHROM, str(rec.POS)] + [(rc[i], dp[i], x) for i, x in enumerate(vaf)]
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
  num_loc = [0] * len(samples)
  sum_depth = [0] * len(samples)
  var_data = []
  for rec in rdr:
    vdat = parse_vcf_record(rec, samples, do_add_normal)
    var_data.append(vdat)
    for i in range(1, len(samples)):
      alt, tot, vaf = vdat[2+i]
      if not tot is None:
        num_loc[i] += 1
        sum_depth[i] += tot
  mean_depth = [round(dp / nloc) if nloc > 0 else 0 for dp, nloc in zip(sum_depth, num_loc)]

  # output header
  hdr = [''] + samples[1:] # first column contains row ids (R style)
  # set up output files
  fn_alt = os.path.join(args.outdir, 'rc_alt.csv')
  f_alt = open(fn_alt, 'wt')
  f_alt.write(','.join(hdr) + '\n')
  fn_tot = os.path.join(args.outdir, 'rc_tot.csv')
  f_tot = open(fn_tot, 'wt')
  f_tot.write(','.join(hdr) + '\n')

  for vdat in var_data:
    chrom, pos = vdat[:2]
    out_dat_alt = ['{}_{}'.format(chrom, pos)]
    out_dat_tot = ['{}_{}'.format(chrom, pos)]
    for i in range(1, len(samples)):
      alt, tot, vaf = vdat[2+i]
      if tot is None: # if variant is not present in this sample, fill in mean depth
        tot = mean_depth[i]
        alt = 0
      out_dat_alt += [str(alt)]
      out_dat_tot += [str(tot)]
    f_alt.write(','.join(out_dat_alt) + '\n')
    f_tot.write(','.join(out_dat_tot) + '\n')

  f_alt.close()
  f_tot.close()

if __name__ == '__main__':
  args = parse_args()
  main(args)
