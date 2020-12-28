#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Parse PhyloWGS output files to extract clusters and tree.
#
# Output:
#  - cluster.csv : id_cluster,id_sample,freq
#  - snv.csv     . chrom_pos,id_cluster,id_snv
#  - tree.csv    . id_tree,from_cluster,to_cluster
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-29
#------------------------------------------------------------------------------

import os, sys
import argparse
import gzip, json, zipfile

def parse_args():
  parser = argparse.ArgumentParser(description='Extract clusters and tree from PhyloWGS output.')
  parser.add_argument('input_snv', type=argparse.FileType(), help='Path to input file with SNV info.')
  parser.add_argument('result_summ', type=argparse.FileType(), help='Path to summary result file (*.summ.json.gz).')
  parser.add_argument('result_mutass', type=argparse.FileType(), help='Path to mutation assignments result file (*.mutass.zip).')
  parser.add_argument('--outdir', required=False, help='Path to output directory (default: same as input).')

  args = parser.parse_args()
  args_dict = vars(args)
  if not args.outdir:
    args_dict['outdir'] = os.path.dirname(args.result_summ.name) or '.'

  return args_dict

def parse_input(fh):
  snvs = {}
  # read header
  hdr = fh.readline()
  # read SNV data
  for line in fh:
    cols = line.strip().split('\t')
    snvs[cols[0]] = cols[1]

  return snvs

def parse_summ(fh):
  trees = []
  d = json.load(gzip.open(fh.name, 'r'))
  # extract log-likelihoods for trees
  llh = [d['trees'][t]['llh'] for t in d['trees']]
  # determine tree with best LLH (choose last if tied)
  max_llh = max(llh)
  idx_max_llh = [i for i, j in enumerate(llh) if j == max_llh][-1]
  # extract cellular prevalence for clusters
  clust_prev = {k: v['cellular_prevalence'] for k, v in d['trees'][str(idx_max_llh)]['populations'].items()}
  # extract clone tree branches for clusters
  clust_children = d['trees'][str(idx_max_llh)]['structure']

  return idx_max_llh, clust_prev, clust_children

def parse_mutass(fh, idx):
  z = zipfile.ZipFile(fh.name)
  d = json.loads(z.read('{}.json'.format(idx))) 
  clust_snv = {k: v['ssms'] for k, v in d['mut_assignments'].items()}
  
  return clust_snv

def main(args):
  print(args)
  variants = parse_input(args['input_snv'])
  idx_best, clust_prev, clust_children = parse_summ(args['result_summ'])
  clust_var = parse_mutass(args['result_mutass'], idx_best)

  print("%d clusters" % len(clust_var))
  print("%d snvs" % len(variants))
  
  # write clusters to file
  fn_clust = os.path.join(args['outdir'], 'inf.clusters.csv')
  with open(fn_clust, 'wt') as f:
    f.write('id_cluster,id_sample,freq\n')
    for id_clust in sorted(clust_prev.keys()):
      freq = clust_prev[id_clust]
      for i, prev in enumerate(freq):
        f.write('{},{},{}\n'.format(id_clust, 'R{}'.format(i), prev))

  # write trees to file
  fn_trees = os.path.join(args['outdir'], 'inf.trees.csv')
  with open(fn_trees, 'wt') as f:
    f.write('from,to\n')
    for id_clust, children in clust_children.items():
      for child in children:
        f.write('{},{}\n'.format(id_clust, child))

  # write variants to file
  fn_snvs = os.path.join(args['outdir'], 'inf.snvs.csv')
  with open(fn_snvs, 'wt') as f:
    f.write('chrom_pos,id_cluster,id_snv\n')
    for id_clust, lst_vars in clust_var.items():
      for id_var in lst_vars:
        chrom_pos = variants[id_var]
        f.write('{},{},{}\n'.format(chrom_pos, id_clust, id_var))

if __name__ == '__main__':
  args = parse_args()
  main(args)
