#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert multi-sample VCF file to PyClone input TSV files.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-20
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import vcf
import math

# output header
#hdr = "id_mut rc_ref rc_alt cn_normal cn_minor cn_major".split()
hdr = "mutation_id ref_counts var_counts normal_cn minor_cn major_cn".split()

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone input files from VCF.')
  parser.add_argument('--vcf', required=True, type=argparse.FileType(), help='Multi-sample VCF file with somatic mutations and read counts. (required FORMAT fields: "AD", "DP", "CNI"')
  parser.add_argument('--normal', required=True, help='Normal sample id. (ignored for output)')
  parser.add_argument('--outdir', required=True, help='Output directory for TSV files.')

  args = parser.parse_args()
  return args

def parse_vcf_record(rec, lbl_samples):
  '''
  Parse a multi-sample VCF record.
  Determine major ALT allele and calculate VAF for each sample.
  
  OUTPUT: 
    list with the following elements:
      CHROM_POS REF ALT (S1:ref, S1:alt, S1:CNI) [(S2:ref, S2:alt, S2:CNI) [...]]
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
  cn = [None] * len(idx_smp)
  for i in idx_smp:
    if i < len(rec.samples):
      cdata = rec.samples[i].data
      ad = getattr(cdata, 'AD') if hasattr(cdata, 'AD') else None
      if ad:
        rc[i] = ad[idx_maj+1]
      dp[i] = getattr(cdata, 'DP')
      cn[i] = getattr(cdata, 'CNI') if hasattr(cdata, 'CNI') else 2
      
  # compile REF/ALT counts and inferred copy number
  alt = [  x if x else 0 for   x in rc]
  ref = [d-a if d else 0 for a,d in zip(alt, dp)]

  desc = rec.ID if rec.ID else '.'
  id_var = '{}_{}'.format(rec.CHROM, str(rec.POS)) 
  out = [id_var, rec.REF, major_allele] + [x if x[2] else (x[0], x[1], 2) for x in zip(ref, alt, cn)]
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
  lst_vars = []
  for rec in rdr:
    var_data = parse_vcf_record(rec, samples)
    lst_vars.append(var_data)

  # write output files (ignore normal sample)
  for i in range(1, len(samples)):
    id_smp = samples[i]
    fn_out = os.path.join(args.outdir, id_smp + '.tsv')
    with open(fn_out, 'wt') as fh_out:
      fh_out.write('\t'.join(hdr) + '\n')
      for var_data in lst_vars:
        id_var = var_data[0]
        ref_alt_cn = var_data[3+i]
        rc_ref = ref_alt_cn[0]
        rc_alt = ref_alt_cn[1]
        cn_tot = ref_alt_cn[2]
        # cn_nrm = var_data[3][2] # normal sample missing CNI frequently
        cn_nrm = 2
        # heuristic to split total CN in minor and major by VAF
        vaf = (rc_alt / (rc_alt + rc_ref)) if (rc_alt + rc_ref) > 0 else 0.0
        cn_min = math.ceil(vaf*cn_tot)
        if rc_ref > 0: # cn_min should be lower than cn_tot if there are REF reads
          cn_min = min(cn_min, cn_tot-1)
        cn_maj = cn_tot - cn_min 
        if True: #rc_alt > 0:
          out = [str(x) for x in [id_var, rc_ref, rc_alt, cn_nrm, cn_min, cn_maj]]
          assert len(out) == len(hdr), 'Number of data fields does not match header.'
          fh_out.write('\t'.join(out) + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args)
