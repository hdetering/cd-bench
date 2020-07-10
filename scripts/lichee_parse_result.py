#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Parse LICHeE output file to extract clusters and tree.
#
# Output:
#  - cluster.csv : id_cluster,presence,mean_freq
#  - snv.csv     . id_snv,chrom,pos,id_cluster
#  - tree.csv    . id_tree,from_cluster,to_cluster
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-05-17
#------------------------------------------------------------------------------

import os, sys
import argparse
import re

def parse_args():
  parser = argparse.ArgumentParser(description='Extract clusters and trees from LICHeE output.')
  parser.add_argument('lichee_result', type=argparse.FileType(), help='Path to LICHeE output file.')
  parser.add_argument('--outdir', required=False, help='Path to output directory (default: same as input).')

  args = parser.parse_args()
  args_dict = vars(args)
  if not args.outdir:
    args_dict['outdir'] = os.path.dirname(args.lichee_result.name) or '.'

  return args_dict

def parse_lichee(fh):
  clusters = {}
  trees = {}
  snvs = {}

  do_parse_clusters = False
  do_parse_trees = False
  do_parse_snvs = False

  id_tree = -1

  for line in fh:
    if line.startswith('Nodes:'):
      print('Parsing clusters...')
      do_parse_clusters = True
    elif line.startswith('****Tree'):
      print('Parsing trees...')
      do_parse_clusters = False
      do_parse_trees = True
      m = re.search('Tree (\d+)', line)
      id_tree = int(m.group(1))
    elif line.startswith('SNV info:'):
      print('Parsing SNVs...')
      do_parse_trees = False
      do_parse_snvs = True
    elif do_parse_clusters:
      cols = line.strip().split('\t')
      if len(cols) > 3:
        id_cluster, pres, freq = cols[:3]
        clusters[id_cluster] = (pres, freq, cols[3:])
    elif do_parse_trees:
      if id_tree not in trees:
        trees[id_tree] = []
      m = re.search('(\d+) -> (\d+)', line)
      if m:
        trees[id_tree].append((m.group(1), m.group(2)))
    elif do_parse_snvs:
      m = re.search('(snv\d+): (\w+) (\d+) ([a-zA-Z0-9]+)$', line)
      if m:
        snvs[m.group(1)] = (m.group(2), m.group(3), m.group(4))
  
  return clusters, trees, snvs 

def main(args):
  print(args)
  clusters, trees, variants = parse_lichee(args['lichee_result'])

  print("%d clusters" % len(clusters))
  print("%d trees" % len(trees))
  print("%d snvs" % len(variants))
  
  # remember which cluster variants are assigned to
  snv_clust = {}

  # write clusters to file
  fn_clust = os.path.join(args['outdir'], 'lichee.clusters.csv')
  with open(fn_clust, 'wt') as f:
    f.write('id_cluster,pres,freq\n')
    for id_clust in sorted(clusters.keys()):
      pres, freq, snvs = clusters[id_clust]
      f.write('{},{},{}\n'.format(id_clust, pres, freq))
      for id_snv in snvs:
        snv_clust[id_snv] = id_clust

  # write trees to file
  fn_trees = os.path.join(args['outdir'], 'lichee.trees.csv')
  with open(fn_trees, 'wt') as f:
    f.write('id_tree,from,to\n')
    for id_tree in sorted(trees.keys()):
      for from_to in trees[id_tree]:
        f.write('{},{},{}\n'.format(id_tree, from_to[0], from_to[1]))

  # write variants to file
  fn_snvs = os.path.join(args['outdir'], 'lichee.snvs.csv')
  with open(fn_snvs, 'wt') as f:
    f.write('id_mut,id_cluster,id_snv,chrom,pos\n')
    for id_snv in sorted(variants.keys()):
      chrom, pos, name = variants[id_snv]
      f.write('{},{},{},{},{}\n'.format(name, snv_clust[id_snv], id_snv, chrom, pos))

if __name__ == '__main__':
  args = parse_args()
  main(args)
