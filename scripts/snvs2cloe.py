#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert TSV file with SNV info to Cloe input files.
#
# INPUT: chrom, pos, smp, ref, alt, rc_tot, rc_alt
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-11-08
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys
import argparse
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Create Cloe input files from CSV.')
  parser.add_argument('--csv', required=True, type=argparse.FileType(), help='CSV file with somatic mutations and read counts.')
  parser.add_argument('--outdir', required=True, help='Output directory.')
  parser.add_argument('--normal', help='Normal sample id. (excluded from output)')

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

  var_smp_rc, set_samples = parse_csv(args.csv)
  samples = [x for x in sorted(set_samples) if x != args.normal]

  # output header
  hdr = [''] + samples # first column contains row ids (R style)
  # set up output files
  fn_alt = os.path.join(args.outdir, 'rc_alt.csv')
  f_alt = open(fn_alt, 'wt')
  f_alt.write(','.join(hdr) + '\n')
  fn_tot = os.path.join(args.outdir, 'rc_tot.csv')
  f_tot = open(fn_tot, 'wt')
  f_tot.write(','.join(hdr) + '\n')

  for id_var in sorted(var_smp_rc.keys()):
    out_dat_alt = [id_var]
    out_dat_tot = [id_var]
    for id_smp in samples:
      smp_rc = var_smp_rc[id_var]
      tot, alt = (0, 0)
      if id_smp in smp_rc:
        tot, alt = smp_rc[id_smp]
      out_dat_alt += [str(alt)]
      out_dat_tot += [str(tot)]
    f_alt.write(','.join(out_dat_alt) + '\n')
    f_tot.write(','.join(out_dat_tot) + '\n')

  f_alt.close()
  f_tot.close()

if __name__ == '__main__':
  args = parse_args()
  main(args)
