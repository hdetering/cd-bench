#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert multi-sample VCF file to SciClone input CSV files.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-11-20
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import vcf
import math

# output header
hdr = "chr pos refCount varCount VAF".split()

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone input files from VCF.')
  parser.add_argument('--csv', required=True, type=argparse.FileType(), help='CSV file with somatic mutations and read counts.')
  parser.add_argument('--normal', help='Normal sample id. (ignored for output)')
  parser.add_argument('--outdir', required=True, help='Output directory for CSV files.')

  args = parser.parse_args()
  return args

def parse_csv(fh):
  '''
  Parse contents of CSV file with the following column specification:
    chrom,pos,id_sample,ref,alt,rc_tot,rc_alt
  Output:
    dict {id_var: {id_sample: (rc_tot, rc_alt)}}
  '''
  samples = set()
  var_smp_rc = {}
  # skip header line
  hdr = fh.readline()
  # parse variant data
  for line in fh:
    cols = line.strip().split(',')
    chrom, pos, smp, ref, alt, rc_tot, rc_alt = cols[:7]
    samples.update((smp,))
    id_var = '{}_{}'.format(chrom, pos)
    if id_var not in var_smp_rc:
      var_smp_rc[id_var] = {}
    assert smp not in var_smp_rc[id_var], 'Variant "{}" found more than once for sample "{}"'.format(id_var, smp)
    var_smp_rc[id_var][smp] = (rc_tot, rc_alt)

  return var_smp_rc, samples

def main(args):
  # read variant data
  var_smp_rc, set_samples = parse_csv(args.csv)
  samples = [x for x in sorted(set_samples) if x != args.normal]

  # write output files (ignore normal sample)
  for i in range(len(samples)):
    id_smp = samples[i]
    fn_out = os.path.join(args.outdir, id_smp + '.vaf.csv')
    with open(fn_out, 'wt') as fh_out:
      fh_out.write(','.join(hdr) + '\n')
      for id_var in sorted(var_smp_rc.keys()):
        chrom, pos = id_var.split('_')
        smp_rc = var_smp_rc[id_var]
        rc_tot, rc_alt = (0, 0)
        if id_smp in smp_rc:
          str_tot, str_alt = smp_rc[id_smp]
          rc_tot = int(str_tot)
          rc_alt = int(str_alt)
        rc_ref = rc_tot - rc_alt
        # heuristic to split total CN in minor and major by VAF
        vaf = (rc_alt / (rc_alt + rc_ref)) if (rc_alt + rc_ref) > 0 else 0.0
        str_vaf = '{:.4f}'.format(vaf)
        out = [str(x) for x in [chrom, pos, rc_ref, rc_alt, str_vaf]]
        assert len(out) == len(hdr), 'Number of data fields does not match header.'
        fh_out.write(','.join(out) + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args)
