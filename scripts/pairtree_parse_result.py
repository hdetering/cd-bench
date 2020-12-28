#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Extract PairTree results and convert to canonical format.
#
# INPUT:
#   .npz archive
# OUTPUT:
#   inf.snvs.csv
#     chrom_pos,id_cluster
#   inf.clusters.csv
#     id_cluster,id_sample,freq
#   inf.trees.csv
#     id_tree,from,to 
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-12-21
#------------------------------------------------------------------------------

import argparse
import numpy as np
import json
import os

fn_out_snvs  = 'inf.snvs.csv'
fn_out_clust = 'inf.clusters.csv'
fn_out_tree  = 'inf.trees.csv'

def parse_args():
  parser = argparse.ArgumentParser(description='Extract PairTree results from .npz archive.')
  parser.add_argument('npz', type=argparse.FileType(), help='PairTree result (.npz archive)')
  parser.add_argument('outdir', help='Output directory.')

  args = parser.parse_args()
  return args

if __name__ == '__main__':
  args = parse_args()

  npz = np.load(args.npz.name)
  # [x for x in npz.keys()]
  # ['accept_rate.json',
  #  'struct',
  #  'count',
  #  'phi',
  #  'llh',
  #  'prob',
  #  'sampnames.json',
  #  'clustrel_posterior_vids.json',
  #  'clustrel_posterior_rels',
  #  'clustrel_evidence_vids.json',
  #  'clustrel_evidence_rels',
  #  'clusters.json',
  #  'garbage.json',
  #  'seed.json']

  lst_clust = json.loads(npz['clusters.json'])
  with open(os.path.join(args.outdir, fn_out_snvs), 'wt') as f:
    f.write('chrom_pos,id_cluster\n')
    for (idx_clust,), lst_muts in np.ndenumerate(lst_clust):
      for id_mut in lst_muts:
        f.write('{},{}\n'.format(id_mut, idx_clust))

  # extract sample names
  samples = json.loads(npz['sampnames.json'])

  # get maximum log likelihood
  llh = list(npz['llh'])
  llh_max = max(llh)
  idx_max = llh.index(llh_max)

  # extract clone prevalence
  phi = npz['phi'][idx_max]
  with open(os.path.join(args.outdir, fn_out_clust), 'wt') as f:
    f.write('id_cluster,id_sample,freq\n')
    for (idx_clust, idx_smp), prev in np.ndenumerate(phi):
      f.write('{},{},{:.4f}\n'.format(idx_clust,samples[idx_smp], prev))

  # extract clone tree
  struct = npz['struct'][idx_max]
  with open(os.path.join(args.outdir, fn_out_tree), 'wt') as f:
    f.write('id_tree,from,to\n')
    for (idx_child,), idx_parent in np.ndenumerate(struct):
      f.write('{},{},{}\n'.format(idx_max, idx_parent, idx_child+1))
