#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Calculate distance measures between two trees.
#
# Input:
#   - TRUE clone tree adjacency list (CSV: from,to)
#   - TRUE variant-to-clone assignment (CSV: id_var,chrom,pos,id_clone)
#   - INFERRED tree (CSV: from,to)
#   - INFERRED variant-to-cluster assignment (CSV: id,chrom,pos,id_var,id_cluster)
# Output:
#   distance value
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-27
#------------------------------------------------------------------------------

from __future__ import division
import argparse
import itertools as it
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Calculate tree distance.')
  # positional arguments
  parser.add_argument('true_tree', type=argparse.FileType('r'), help='TRUE clone tree adjacency list (CSV: from,to).')
  parser.add_argument('true_snvs', type=argparse.FileType('r'), help='TRUE variant-to-clone mapping (CSV: id_var,chrom,pos,id_clone).')
  parser.add_argument('inf_tree', type=argparse.FileType('r'), help='INFERRED tree (CSV: from,to).')
  parser.add_argument('inf_snvs', type=argparse.FileType('r'), help='INFERRED variant-to-cluster mapping (CSV: id,chrom,pos,id_var,id_cluster).')
  # optional arguments
  grp_metric = parser.add_mutually_exclusive_group(required=True)
  grp_metric.add_argument('--CASet', action='store_true')
  grp_metric.add_argument('--DISC', action='store_true')
  grp_variant = parser.add_mutually_exclusive_group()
  grp_variant.add_argument('--intersect', action='store_true')
  grp_variant.add_argument('--union', action='store_true')
  
  args = parser.parse_args()
  return args

def Jaccard_dist(A, B):
  '''
  Return the Jaccard distance between sets A and B.
  '''
  assert type(A) is set and type(B) is set, "A and B must be sets."
  jacc_dist = (len(A|B) - len(A&B)) / len(A|B) if len(A|B) > 0 else 0.0
  return jacc_dist

def A(id_var, tree, clusters, var_clust):
  '''
  Return ancestral set for a given variant.

  The ancestral set A(i) includes all variants from the root to the cluster
  which contains i. I.e., it is the set of mutations that make up the clone
  in which i first occurred.
  '''
  # cluster which variant belongs to
  clust = var_clust[id_var]
  # determine all ancestral clusters
  ancestors = [clust]
  clust_node = None
  for node in tree.iter_descendants():
    if node.name == clust:
      clust_node = node
      break
  if clust_node is None:
    import pdb; pdb.set_trace()
  assert clust_node, "Cluster not found in tree."
  ancestors += [n.name for n in clust_node.iter_ancestors()]
  # collect variants in ancestral clusters
  anc_vars = set()
  for c in ancestors:
    if c in clusters: # empty clusters are not in dict
      anc_vars.update(clusters[c])
  
  return anc_vars

def CASet(variants, tree1, tree2, clust1, clust2, snv_clust1, snv_clust2):
  '''
  Calculate Common Ancestor Set (CASet) distance between two trees.
  '''
  n_comp = 0 # number of comparisons
  cum_dist = 0.0 # cumulative distance
  for i, j in it.combinations(variants, 2):
    n_comp += 1
    # calculate ancestral sets for mutations i and j in tree1
    A1_i = A(i, tree1, clust1, snv_clust1) if i in snv_clust1 else set()
    A1_j = A(j, tree1, clust1, snv_clust1) if j in snv_clust1 else set()
    # calculate common ancestor set of mutations i and j in tree1
    C1 = A1_i & A1_j
    # calculate ancestral sets for mutations i and j in tree2
    A2_i = A(i, tree2, clust2, snv_clust2) if i in snv_clust2 else set()
    A2_j = A(j, tree2, clust2, snv_clust2) if j in snv_clust2 else set()
    # calculate common ancestor set of mutations i and j in tree2
    C2 = A2_i & A2_j

    jacc_dist = Jaccard_dist(C1, C2)
    cum_dist += jacc_dist

  return cum_dist / n_comp

def adjacency_list_to_tree(adj_list):
  '''Initialize ete3.Tree using a set of edges.'''
  nodes = {}
  #import pdb; pdb.set_trace()
  for lbl_parent, lbl_child in adj_list:
    p = nodes[lbl_parent] if lbl_parent in nodes else Tree(name=lbl_parent)
    c = nodes[lbl_child] if lbl_child in nodes else Tree(name=lbl_child)
    p.add_child(c)
    nodes[lbl_parent] = p
    nodes[lbl_child] = c
  
  roots = [n for lbl, n in nodes.items() if n.is_root()]
  assert len(roots) == 1, "More than one root found."
  return roots[0]

def main(args):
  true_tree_edges = []
  true_snv_clust = {}
  true_clusters = {}
  inf_tree_edges = []
  inf_snv_clust = {}
  inf_clusters = {}

  # read input files
  #-----------------------------------------------------------------------------
  for line in args.true_tree:
    row = line.strip().split(',')
    true_tree_edges.append((row[0],row[1]))
  line = next(args.true_snvs) # skip header row
  for line in args.true_snvs:
    row = line.strip().split(',')
    id_snv = row[0]
    id_clust = row[3]
    true_snv_clust[id_snv] = id_clust
    if id_clust in true_clusters:
      true_clusters[id_clust].append(id_snv)
    else:
      true_clusters[id_clust] = [id_snv]
  for line in args.inf_tree:
    row = line.strip().split(',')
    inf_tree_edges.append((row[0],row[1]))
  line = next(args.inf_snvs) # skip header row
  for line in args.inf_snvs:
    row = line.strip().split(',')
    id_snv = row[3]
    id_clust = row[4]
    inf_snv_clust[id_snv] = id_clust
    if id_clust in inf_clusters:
      inf_clusters[id_clust].append(id_snv)
    else:
      inf_clusters[id_clust] = [id_snv]
  # convert edge lists to trees
  true_tree = adjacency_list_to_tree(true_tree_edges)  
  inf_tree = adjacency_list_to_tree(inf_tree_edges)  

  # sanity checks
  #-----------------------------------------------------------------------------
  # check that ISA holds
  true_vars_lst = true_snv_clust.keys()
  true_vars_set = set(true_vars_lst)
  assert len(true_vars_set) == len(true_vars_lst), "TRUE tree violates ISA."
  inf_vars_lst = inf_snv_clust.keys()
  inf_vars_set = set(inf_vars_lst)
  assert len(inf_vars_set) == len(inf_vars_lst), "INFERRED tree violates ISA."
  if not args.intersect and not args.union:
    # trees must contain identical variant sets
    assert set(true_vars) == set(inf_vars), "Trees contain different variant sets. (use --intersect/--union)"
  
  # calculate distance
  #-----------------------------------------------------------------------------
  tree_dist = None
  if args.CASet:
    # compile variant set
    snvs = true_vars_set
    if args.intersect:
      snvs = true_vars_set & inf_vars_set
    elif args.union:
      snvs = true_vars_set | inf_vars_set
    # perform distance calculation
    tree_dist = CASet(snvs, true_tree, inf_tree, true_clusters, inf_clusters, true_snv_clust, inf_snv_clust)

  print('Tree distance: {:.4f}'.format(tree_dist))

if __name__ == '__main__':
  args = parse_args()
  main(args)
