#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert SNV infos to PairTree input files.
#
# INPUT:
#   CSV format: chrom,pos,id_sample,ref,alt,rc_tot,rc_alt
# OUTPUT:
#   - TSV format: id gene a d mu_r mu_v
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-12-19
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys
import argparse
import json

def parse_args():
  parser = argparse.ArgumentParser(description='Create PairTree input files from SNV file.')
  parser.add_argument('--csv', required=True, type=argparse.FileType(), help='CSV file with somatic mutations and read counts.')
  parser.add_argument('--outdir', required=True, help='Directory to write output files to (must exist).')
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

  fh_out = open(os.path.join(args.outdir, 'input.snvs.tsv'), 'wt')
  # output header
  hdr = "id name var_reads total_reads var_read_prob".split()
  fh_out.write('\t'.join(hdr) + '\n')

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
    prb = ','.join(['0.5']*len(samples))
    # import pdb; pdb.set_trace()
    #alt = ','.join([str(a) if a else '0' for s in samples for t,a in smp_rc[s]])
    #tot = ','.join([str(t) if t else '0' for s in samples for t,a in smp_rc[s]])
    out_data = [id_mut, id_var, alt, tot, prb]
    fh_out.write('\t'.join(out_data) + '\n')
    idx_var += 1

  # write params file
  with open(os.path.join(args.outdir, 'params.pre.json'), 'wt') as f:
    f.write(json.dumps({"samples": samples}) + '\n')


if __name__ == '__main__':
  args = parse_args()
  main(args)
