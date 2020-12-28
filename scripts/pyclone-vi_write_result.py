#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert PyClone-VI H5 output file to two TSV files.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-10-14
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import pyclone_vi.post_process

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone YAML input file from VCF and copy-number BED.')
  parser.add_argument('--h5', required=True, type=argparse.FileType(), help='Output of pyclone-vi fit (H5 file).')
  parser.add_argument('--outdir', required=True, help='Output directory for TSV files.')

  args = parser.parse_args()
  return args


def main(fh_input, outdir):
  '''Load data from input and split into two files.'''
  df_clust = pyclone_vi.post_process.load_cluster_df(fh_input.name)
  df_clust = df_clust.rename(columns={'sample_id': 'id_sample', 'cluster_id': 'id_cluster', 'cellular_prevalence': 'freq'})
  df_clust = df_clust[['id_cluster', 'id_sample', 'freq']]
  df_loci = pyclone_vi.post_process.load_loci_df(fh_input.name)
  df_loci = df_loci.rename(columns={'mutation_id': 'chrom_pos', 'cluster_id': 'id_cluster'})
  df_loci = df_loci[['chrom_pos', 'id_cluster']]

  fn_out_clust = os.path.join(outdir, 'inf.clusters.csv')
  fn_out_loci  = os.path.join(outdir, 'inf.snvs.csv')

  df_clust.to_csv(fn_out_clust, float_format='%.4f', index=False, sep=',')
  df_loci.to_csv(fn_out_loci, float_format='%.4f', index=False, sep=',')

if __name__ == '__main__':
  args = parse_args()
  main(args.h5, args.outdir)
