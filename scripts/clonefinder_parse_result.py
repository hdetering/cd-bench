#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Parse CloneFinder output file to extract clusters and prevalences.
#
# Output:
#  - inf.clusters.csv : id_sample,id_cluster,freq
#  - inf.snvs.csv     : id_snv,id_cluster
#  - inf.trees.csv    : id_tree,from,to
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-11-12
#------------------------------------------------------------------------------

import os, sys
import argparse
#import re
from bitarray import bitarray
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Extract clusters and tree from CloneFinder output.')
  parser.add_argument('input', type=argparse.FileType(), help='Path to CloneFinder input file (.tsv).')
  parser.add_argument('tree', type=argparse.FileType(), help='Path to CloneFinder tree file (.nwk).')
  parser.add_argument('genotypes', type=argparse.FileType(), help='Path to CloneFinder genotypes file (.meg).')
  parser.add_argument('prevalences', type=argparse.FileType(), help='Path to CloneFinder prevalences file (.txt).')
  parser.add_argument('--outdir', required=False, help='Path to output directory (default: same as input).')

  args = parser.parse_args()
  args_dict = vars(args)
  if not args.outdir:
    args_dict['outdir'] = os.path.dirname(args.tree.name) or '.'

  return args_dict

def parse_input_vars(fh):
  loci = []
  # read header line
  line = fh.readline()
  # parse variant loci
  for line in fh:
    parts = line.strip().split('\t')
    loci.append(parts[0])
  
  return loci

def parse_tree(fh):
  '''Parse CloneFinder tree file.'''
  t = Tree(fh.name, format=1)
  # make sure all nodes are labelled
  i = 0
  for node in t.traverse('preorder'):
    if not node.name:
      node.name = 'node_{}'.format(i)
      i += 1

  return t

def get_genotype(node, node_gt):
  if node.name in node_gt:
    gt = node_gt[node.name]
  else: # reconstruct GT as intersection of children's GTs
    lst_child_gt = []
    for child in node.children:
      lst_child_gt.append(get_genotype(child, node_gt))
    gt = lst_child_gt[0] # if this fails, we arrived at a tip without GT -> baaad!
    if len(lst_child_gt) > 1: # if there is just one child, parent is assumed to share same GT
      for child_gt in lst_child_gt[1:]:
        gt = gt & child_gt
    # subtract parent muts from child GT (record only first occurrence of muts)
    for child in node.children:
      child_gt = node_gt[child.name]
      node_gt[child.name] = child_gt &~ gt
    node_gt[node.name] = gt

  return gt

def reduce_genotype(node, node_gt):
  '''Recursively remove mutations of parent from child nodes.'''
  gt_parent = node_gt[node.name]
  for child in node.children:
    gt_child = node_gt[child.name]
    # update child genotype (gt_child AND (NOT gt_parent))
    node_gt[child.name] = gt_child &~ gt_parent
    reduce_genotype(child, node_gt)

def parse_genotypes(fh, tree):
  '''
  Parse CloneFinder genotypes file.
  
  Genotypes are converted, so that mutations are assigned to the tree node
  in which they first ocurred.
  '''
  # number of lines to discard from top of file
  skip = 4
  # skip n lines
  for i in range(skip):
    line = fh.readline()

  node_gt = {} # store genotype (bitarray) for each node
  header = ''
  # parse header/genotype records
  for line in fh:
    if line.startswith('#'):
      header = line.strip()
      continue
    if header == '#hg19':
      # fix header (sample id different from tree node lbl...)
      header = '#Normal'
    # extract node id from header
    id_node = header[1:]
    # parse genotype
    gt_str = line.strip() # expected: string of 'A's (absent) and 'T's (present)
    gt = bitarray([1 if c == 'T' else 0 for c in gt_str])
    node_gt[id_node] = gt
  
  # reconstruct genotypes of internal tree nodes
  root_gt = get_genotype(tree, node_gt)
  assert sum(root_gt) == 0, 'Root genotype should be all REF...'

  # convert genotypes into mutation clusters (only first occurrence of mut noted)
  #import pdb; pdb.set_trace()
  #reduce_genotype(tree, node_gt)

  return node_gt

def parse_prevalence(fh):
  '''Read prevalence matrix and convert to dict.'''
  # read header line
  header = fh.readline().strip().split('\t')
  # set up output dict
  clust_smp_prev = {x: {} for x in header[1:]}
  # read data
  for line in fh:
    cols = line.strip().split()
    id_smp = cols[0]
    for i in range(1, len(cols)):
      id_clust = header[i]
      clust_smp_prev[id_clust][id_smp] = float(cols[i])

  return clust_smp_prev

def main(args):
  #print(args)
  
  # read variant loci from input file
  lst_var = parse_input_vars(args['input'])
  # read clone tree
  tree = parse_tree(args['tree'])
  # reconstruct genotypes
  node_gt = parse_genotypes(args['genotypes'], tree)
  # read cluster prevalences
  clust_smp_prev = parse_prevalence(args['prevalences'])
  
  print("%d clusters" % len(clust_smp_prev))
  #print("%d trees" % len(trees))
  print("%d variant loci" % len(lst_var))
  
  # remember which cluster variants are assigned to
  snv_clust = {}

  # write clusters to file
  fn_clust = os.path.join(args['outdir'], 'inf.clusters.csv')
  with open(fn_clust, 'wt') as f:
    f.write('id_cluster,id_sample,freq\n')
    for id_clust in sorted(clust_smp_prev.keys()):
      for id_smp in sorted(clust_smp_prev[id_clust].keys()):
        freq = clust_smp_prev[id_clust][id_smp]
        f.write('{},{},{:.4f}\n'.format(id_clust, id_smp, freq))

  # write trees to file
  fn_trees = os.path.join(args['outdir'], 'inf.trees.csv')
  with open(fn_trees, 'wt') as f:
    f.write('id_tree,from,to\n')
    for node in tree.traverse('preorder'):
      for child in node.children:
        f.write('0,{},{}\n'.format(node.name, child.name))

  # write variants to file
  fn_snvs = os.path.join(args['outdir'], 'inf.snvs.csv')
  with open(fn_snvs, 'wt') as f:
    f.write('chrom_pos,id_cluster\n')
    for id_cluster in sorted(node_gt.keys()):
      for id_mut, present in zip(lst_var, node_gt[id_cluster]):
        if present:
          f.write('{},{}\n'.format(id_mut, id_cluster))

if __name__ == '__main__':
  args = parse_args()
  main(args)
