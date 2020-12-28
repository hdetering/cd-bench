#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert read count info from multisample VCF file to CSV format.
#
# INPUT:
#   VCF with FORMAT fields: AD, DP, CNIA, CNIB, CN1, CN2
# OUTPUT:
#   chrom,pos,id_sample,ref,alt,rc_tot,rc_alt,cn_maj_inf,cn_min_inf,cn_a_true,cn_b_true
#
# For samples which do not carry the mutation, rc_tot is reported as
# the mean DP of that sample.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-09-08
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys
import argparse
import vcf

out_hdr = 'chrom pos id_sample ref alt rc_tot rc_alt cn_maj_inf cn_min_inf cn_a_true cn_b_true'.split()

def parse_args():
  parser = argparse.ArgumentParser(description='Create Cloe input files from VCF.')
  parser.add_argument('vcf', type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP"')
  #parser.add_argument('--outdir', required=True, help='Output directory.')
  #parser.add_argument('--normal', help='Normal sample id. (default: assume VAF=0.0 for Normal)')

  args = parser.parse_args()
  return args

def parse_vcf_record(rec, lbl_samples):
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
  cnia = [None] * len(idx_smp)
  cnib = [None] * len(idx_smp)
  cn1 = [None] * len(idx_smp)
  cn2 = [None] * len(idx_smp)
  for i in idx_smp:
    if i < len(rec.samples):
      cdata = rec.samples[i].data
      ad = getattr(cdata, 'AD') if hasattr(cdata, 'AD') else None
      if ad:
        rc[i] = ad[idx_maj+1]
      dp[i] = getattr(cdata, 'DP')
      cnia[i] = getattr(cdata, 'CNIA') if hasattr(cdata, 'CNIA') else None
      cnib[i] = getattr(cdata, 'CNIB') if hasattr(cdata, 'CNIB') else None
      cn1[i] = getattr(cdata, 'CN1') if hasattr(cdata, 'CN2') else None
      cn2[i] = getattr(cdata, 'CN2') if hasattr(cdata, 'CN2') else None
      
  # calculate VAFs
  vaf = [float(a)/float(d) if (a and d) else 0.0 for a, d in zip(rc, dp)]

  out  = [rec.CHROM, str(rec.POS), rec.REF, str(major_allele)] 
  out += [(dp[i], rc[i], x, cnia[i], cnib[i], cn1[i], cn2[i]) for i, x in enumerate(vaf)]
  return out

def main(args):
  # get samples from input VCF
  rdr = vcf.Reader(filename=args.vcf.name)
  samples = rdr.samples

  # parse variants
  num_loc = [0] * len(samples)
  sum_depth = [0] * len(samples)
  var_data = []
  for rec in rdr:
    vdat = parse_vcf_record(rec, samples)
    var_data.append(vdat)
    for i in range(len(samples)):
      dp, rc_alt, vaf, cnia, cnib, cn1, cn2 = vdat[4+i]
      if not dp is None:
        num_loc[i] += 1
        sum_depth[i] += dp
  mean_depth = [round(dp / nloc) if nloc > 0 else 0 for dp, nloc in zip(sum_depth, num_loc)]

  # write header
  print(','.join(out_hdr))

  for vdat in var_data:
    chrom, pos, ref, alt,  = vdat[:4]
    out_dat_alt = ['{}_{}'.format(chrom, pos)]
    out_dat_tot = ['{}_{}'.format(chrom, pos)]
    for i, id_smp in enumerate(samples):
      dp, rc_alt, vaf, cnia, cnib, cn1, cn2 = vdat[4+i]
      if dp is None: # if variant is not present in this sample, fill in mean depth
        dp = mean_depth[i]
        rc_alt = 0
      out_dat = [chrom, pos, id_smp, ref, alt, str(dp), str(rc_alt), str(cnia), str(cnib), str(cn1), str(cn2)]
      print(','.join(out_dat))

if __name__ == '__main__':
  args = parse_args()
  main(args)
