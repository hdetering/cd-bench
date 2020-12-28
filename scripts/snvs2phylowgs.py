#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert SNV infos to PhyloWGS input files.
#
# INPUT:
#   CSV format: chrom,pos,id_sample,ref,alt,rc_tot,rc_alt
# OUTPUT:
#   - TSV format: id gene a d mu_r mu_v
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-11-18
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys
import argparse

def parse_args():
  parser = argparse.ArgumentParser(description='Create Canopy input files from SNV file.')
  parser.add_argument('--csv', required=True, type=argparse.FileType(), help='CSV file with somatic mutations and read counts.')
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
  hdr = "id gene a d mu_r mu_v".split()
  print('\t'.join(hdr))

  # write output vars
  idx_var = 0
  for id_var in sorted(var_smp_rc.keys()):
    id_mut = 's{}'.format(idx_var)
    smp_rc = var_smp_rc[id_var]
    lst_alt = []
    lst_tot = []
    for s in samples:
      t, a = smp_rc[s]
      lst_alt.append(a)
      lst_tot.append(t)
    alt = ','.join(lst_alt)
    tot = ','.join(lst_tot)
    # import pdb; pdb.set_trace()
    #alt = ','.join([str(a) if a else '0' for s in samples for t,a in smp_rc[s]])
    #tot = ','.join([str(t) if t else '0' for s in samples for t,a in smp_rc[s]])
    out_data = [id_mut, id_var, alt, tot, '0.999', '0.499']
    print('\t'.join(out_data))
    idx_var += 1


if __name__ == '__main__':
  args = parse_args()
  main(args)
