#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert multi-sample VCF to PhyloWGS input file.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-24
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
  parser.add_argument('--normal', help='Normal sample id. (default: assume VAF=0.0 for Normal)')

  args = parser.parse_args()
  return args

def parse_vcf_record(rec, lbl_samples, add_normal):
  '''
  Parse a multi-sample VCF record.
  Determine major ALT allele and calculate VAF for each sample.
  
  OUTPUT: 
    list with the following columns:
      #chr position (alt_count1, tot_count1), ...
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
      
  out = [rec.CHROM, str(rec.POS)] + [(rc[i], dp[i]) for i in idx_smp]
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
  var_data = []
  smp_dp   = [0] * len(samples) # store cumulative depth per sample
  smp_nvar = [0] * len(samples) # store number of variants per sample
  for rec in rdr:
    vdat = parse_vcf_record(rec, samples, do_add_normal)
    var_data.append(vdat)
    for i in range(len(samples)):
      dp = vdat[2+i][1]
      if not dp is None:
        smp_dp[i] += dp
        smp_nvar[i] += 1
  # calculate mean depth per sample
  smp_mean_dp = [round(dp/nvar) if nvar>0 else 0 for dp, nvar in zip(smp_dp, smp_nvar)]

  # construct header
  fh_out = args.out
  hdr = "id gene a d mu_r mu_v".split()
  fh_out.write('\t'.join(hdr) + '\n')

  # write output vars
  idx_var = 0
  for vdat in var_data:
    if True: #len(var_line) == len(hdr):
      id_mut = 's{}'.format(idx_var)
      id_gene = '{}_{}'.format(vdat[0], vdat[1])
      alt = ','.join([str(a) if a else '0' for a, d in vdat[2:]][1:])
      tot = ','.join([str(a_d[1]) if a_d[1] else str(smp_mean_dp[i]) for i,a_d in enumerate(vdat[2:])][1:])
      out_data = [id_mut, id_gene, alt, tot, '0.999', '0.499']
      fh_out.write('\t'.join(out_data) + '\n')
      idx_var += 1
  

if __name__ == '__main__':
  args = parse_args()
  main(args)
