#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Reduce clone genotypes to cluster genotypes using clone tree.
#
# Background: Some CD tools (e.g. Cloe) report clone genotypes, rather than
# cluster genotypes. That is, ancestral mutations are assigned multiple times,
# to each clone that harbors them. The canonical data type used to compare
# tools used cluster genotypes, i.e. each mutation is assigned to the cluster
# node in the tree where the mutation first arose.
#
# Input:
#  - inf.snvs.csv : chrom_pos,id_clone
#  - inf.tree.csv : id_tree,from,to
#
# Output:
#  id_snv,id_cluster
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-12-22
#------------------------------------------------------------------------------

import os, sys
import argparse
#import re
from bitarray import bitarray
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Extract clusters and tree from CloneFinder output.')
  parser.add_argument('snvs', type=argparse.FileType(), help='Inferred SNV-to-clone assignments.')
  parser.add_argument('tree', type=argparse.FileType(), help='Inferred clone tree.')
  parser.add_argument('--outdir', required=False, help='Path to output directory (default: same as input).')

  args = parser.parse_args()
  args_dict = vars(args)
  if not args.outdir:
    args_dict['outdir'] = os.path.dirname(args.tree.name) or '.'

  return args_dict

def parse_snvs(fh):
  '''Determine clone genotypes from SNV-to-clone assignments.'''
  set_muts = set()
  clone_muts = {}
  # read header line
  line = fh.readline()
  # parse variant loci
  for line in fh:
    row = line.strip().split(',')
    id_mut   = row[0]
    id_clone = row[1]
    set_muts.add(id_mut)
    if id_clone not in clone_muts:
      clone_muts[id_clone] = []
    clone_muts[id_clone].append(id_mut)
  # infer clone genotypes
  lst_muts = list(set_muts)
  clone_gt = {}
  for id_clone, muts in clone_muts.items():
    clone_gt[id_clone] = bitarray([1 if m in muts else 0 for m in lst_muts])

  return lst_muts, clone_gt

def parse_tree(fh):
  '''Initialize ete3.Tree using a set of edges.'''
  nodes = {}
  # parse header line
  hdr = fh.readline().strip().split(',')
  idx_parent = hdr.index('from')
  idx_child  = hdr.index('to') 
  # parse edges
  for line in fh:
    row = line.strip().split(',')
    lbl_parent = row[idx_parent]
    lbl_child  = row[idx_child]
    p = nodes[lbl_parent] if lbl_parent in nodes else Tree(name=lbl_parent)
    c = nodes[lbl_child] if lbl_child in nodes else Tree(name=lbl_child)
    p.add_child(c)
    nodes[lbl_parent] = p
    nodes[lbl_child] = c

  roots = [n for lbl, n in nodes.items() if n.is_root()]
  lbl_nodes = [lbl for lbl, n in nodes.items()]
  assert len(roots) == 1, "More than one root found."
  
  return roots[0], lbl_nodes

def reduce_genotype(node, node_gt):
  '''Recursively remove mutations of parent from child nodes.'''
  gt_parent = node_gt[node.name]
  # call recursively for children
  for child in node.children:
    reduce_genotype(child, node_gt)
    gt_child = node_gt[child.name]
    # update child genotype (gt_child AND (NOT gt_parent))
    node_gt[child.name] = gt_child &~ gt_parent

def main(args):
  # read mutations from input file
  set_mut, clone_gt = parse_snvs(args['snvs'])
  # read clone tree
  tree, lbl_nodes = parse_tree(args['tree'])
  # make sure every node has a genotype assigned (root may be missing)
  for lbl in lbl_nodes:
    if lbl not in clone_gt:
      clone_gt[lbl] = bitarray([0]*len(set_mut))

  #print("%d clusters" % len(clone_gt))
  #print("%d variant loci" % len(set_mut))
  
  # reduce genotypes from clone to cluster level
  reduce_genotype(tree, clone_gt)

  # write reduced mutation-to-cluster assignment
  #fn_snvs = os.path.join(args['outdir'], 'inf.snvs.csv')
  #with open(fn_snvs, 'wt') as f:
  with sys.stdout as f:
    f.write('chrom_pos,id_cluster\n')
    for id_cluster in sorted(clone_gt.keys()):
      for id_mut, present in zip(set_mut, clone_gt[id_cluster]):
        if present:
          f.write('{},{}\n'.format(id_mut, id_cluster))

if __name__ == '__main__':
  args = parse_args()
  main(args)
