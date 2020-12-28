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
#   distance value(s)
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-09-06
#------------------------------------------------------------------------------

from __future__ import division
import sys
import argparse
import itertools as it
import math
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Calculate tree distance.')
  # positional arguments
  parser.add_argument('true_tree', type=argparse.FileType('r'), help='TRUE clone tree adjacency list (CSV: from,to).')
  parser.add_argument('true_snvs', type=argparse.FileType('r'), help='TRUE variant-to-clone mapping (CSV: id_var,chrom,pos,id_clone).')
  parser.add_argument('inf_tree', type=argparse.FileType('r'), help='INFERRED tree (CSV: from,to).')
  parser.add_argument('inf_snvs', type=argparse.FileType('r'), help='INFERRED variant-to-cluster mapping (CSV: id,chrom,pos,id_var,id_cluster).')
  # optional arguments
  #grp_metric = parser.add_mutually_exclusive_group(required=True)
  parser.add_argument('--CASet_isect', action='store_true')
  parser.add_argument('--CASet_union', action='store_true')
  parser.add_argument('--DISC_isect', action='store_true')
  parser.add_argument('--DISC_union', action='store_true')
  parser.add_argument('--collapse', action='store_true', help='Whether to collapse co-clustered mutations (optimization).')
  parser.add_argument('--ignore-homoplasy', action='store_true', help='Whether to remove homoplasious mutations.')
  
  args = parser.parse_args()
  return args

def choose(n, k):
  '''
  Calculate binomial coefficient.
  '''
  f = math.factorial
  # use integer division to avoid overflows
  return f(n) // f(k) // f(n-k)

def Jaccard_dist(A, B):
  '''
  Return the Jaccard distance between sets A and B.
  '''
  assert type(A) is set and type(B) is set, "A and B must be sets."
  jacc_dist = (len(A|B) - len(A&B)) / len(A|B) if len(A|B) > 0 else 0.0
  return jacc_dist

def Jaccard_dist_weighted(A, B):
  '''
  Return the Jaccard distance between *weighted* sets A and B.
  
  Each set element is expected to be a tuple of the form (id, weight).
  Cardinalities are calculated as the sum of element weights.
  '''
  assert type(A) is set and type(B) is set, "A and B must be sets."
  card_union = sum([v for k, v in A|B])
  card_isect = sum([v for k, v in A&B])
  jacc_dist = (card_union - card_isect) / card_union if card_union > 0 else 0.0
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
  for node in tree.traverse('preorder'):
    if node.name == clust:
      clust_node = node
      break
  if clust_node is None:
    #import pdb; pdb.set_trace()
    pass
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
  cum_dist = 0.0 # cumulative distance
  for i, j in it.combinations(variants, 2):
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

    jacc_dist = Jaccard_dist_weighted(C1, C2)
    cum_dist += jacc_dist * (i[1]*j[1])

  # in case of collapsed subclusters, add self-comparisons
  for i in variants:
    if i[1] == 1: # not a collapsed subcluster
      continue
    # calculate ancestral sets for mutations within i in tree1
    A1_i = A(i, tree1, clust1, snv_clust1) if i in snv_clust1 else set()
    # calculate ancestral sets for mutations within i in tree2
    A2_i = A(i, tree2, clust2, snv_clust2) if i in snv_clust2 else set()
    jacc_dist = Jaccard_dist_weighted(A1_i, A2_i)
    cum_dist += jacc_dist * choose(i[1], 2)

  n_vars = sum([v for k, v in variants]) # global number of variants
  n_comp = choose(n_vars, 2) # number of comparisons
  return cum_dist / n_comp

def DISC(variants, tree1, tree2, clust1, clust2, snv_clust1, snv_clust2):
  '''
  Calculate Distinctly Inherited Set Comparison (DISC) distance between two trees.
  '''
  cum_dist = 0.0 # cumulative distance
  for i in variants:
    for j in variants:
      if i == j:
        continue 
      # calculate ancestral sets for mutations i and j in tree1
      A1_i = A(i, tree1, clust1, snv_clust1) if i in snv_clust1 else set()
      A1_j = A(j, tree1, clust1, snv_clust1) if j in snv_clust1 else set()
      # calculate common ancestor set of mutations i and j in tree1
      D1 = A1_i.difference(A1_j)
      # calculate ancestral sets for mutations i and j in tree2
      A2_i = A(i, tree2, clust2, snv_clust2) if i in snv_clust2 else set()
      A2_j = A(j, tree2, clust2, snv_clust2) if j in snv_clust2 else set()
      # calculate common ancestor set of mutations i and j in tree2
      D2 = A2_i.difference(A2_j)
      
      jacc_dist = Jaccard_dist_weighted(D1, D2)
      cum_dist += jacc_dist * (i[1]*j[1])

  n_vars = sum([v for k, v in variants]) # global number of variants
  n_comp = n_vars * (n_vars-1) # number of comparisons
  return cum_dist / n_comp

def parse_tree_csv(fh):
  '''
  Parse a CSV file containing a set of edges.
  
  Perform the following steps:
    - check if header line is present (expected col names: "from", "to")
      --> extract columns "from", "to"
    - if no header is present, check number of columns
       2 cols --> assume cols "from", "to"
      >2 cols --> assume cols "id_tree", "from", "to"

  Returns: list of tuples [(from_1, to_1), (from_2, to_2), ...]  
  '''
  char_sep = ','
  idx_from = 0
  idx_to = 1
  header = False

  line = fh.readline()
  # check column format
  cols = line.strip().split(char_sep)
  assert len(cols) > 1, "Failed to load adjacency list from %s" % fh.name
  # check for header
  try:
    idx_from = cols.index('from')
    idx_to = cols.index('to')
    header = True
  except: # maybe there's no header, but that's ok
    pass
  if not header:
    # reset file position
    fh.seek(0)
    if len(cols) > 2:
      idx_from = 1
      idx_to = 2
  
  # parse file, populate adjacency list
  from_to = []
  for line in fh:
    cols = line.strip().split(char_sep)
    from_to.append((cols[idx_from], cols[idx_to]))

  return from_to 

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
  inf_snv_cnt = {} # count occurrences of mutations (for homoplasy handling)
  inf_clusters = {}

  # read input files
  #-----------------------------------------------------------------------------
  #for line in args.true_tree:
  #  row = line.strip().split(',')
  #  true_tree_edges.append((row[0],row[1]))
  true_tree_edges = parse_tree_csv(args.true_tree)
  line = next(args.true_snvs) # skip header row
  for line in args.true_snvs:
    row = line.strip().split(',')
    id_snv = row[0]
    id_clust = row[1]
    v = (id_snv, 1) # store cardinality as property of variant
    true_snv_clust[v] = id_clust # store cluster
    if id_clust in true_clusters:
      true_clusters[id_clust].append(v)
    else:
      true_clusters[id_clust] = [v]
  #for line in args.inf_tree:
  #  row = line.strip().split(',')
  #  inf_tree_edges.append((row[0],row[1]))
  inf_tree_edges = parse_tree_csv(args.inf_tree)
  line = next(args.inf_snvs) # skip header row
  for line in args.inf_snvs:
    row = line.strip().split(',')
    id_snv = row[0]
    id_clust = row[1]
    v = (id_snv, 1) # store cardinality as property of variant
    if not args.ignore_homoplasy:
      assert v not in inf_snv_clust, "INFERRED tree violates ISA."
    else:
      if v not in inf_snv_cnt:
        inf_snv_cnt[v] = 1
      else:
        inf_snv_cnt[v] += 1
    inf_snv_clust[v] = id_clust # store cluster
    if id_clust in inf_clusters:
      inf_clusters[id_clust].append(v)
    else:
      inf_clusters[id_clust] = [v]

  # remove homoplasy mutations
  if args.ignore_homoplasy:
    muts_to_remove = [k for k,v in inf_snv_cnt.items() if v>1]
    for var in muts_to_remove:
      del inf_snv_clust[var]
    for id_clust in inf_clusters.keys():
      inf_clusters[id_clust] = [x for x in inf_clusters[id_clust] if x not in muts_to_remove]

  # convert edge lists to trees
  true_tree = adjacency_list_to_tree(true_tree_edges)  
  inf_tree = adjacency_list_to_tree(inf_tree_edges)  

  # sanity checks
  #-----------------------------------------------------------------------------
  # check that ISA holds
  true_vars_lst = true_snv_clust.keys()
  true_vars_set = set(true_vars_lst)
  #assert len(true_vars_set) == len(true_vars_lst), "TRUE tree violates ISA."
  inf_vars_lst = inf_snv_clust.keys()
  inf_vars_set = set(inf_vars_lst)
  #assert len(inf_vars_set) == len(inf_vars_lst), "INFERRED tree violates ISA."
  
  print('TRUE vars: {}'.format(len(true_vars_lst)), file=sys.stderr)
  print('INF vars:  {}'.format(len(inf_vars_lst)), file=sys.stderr)

  # collapse co-clustered variants if requested
  #-----------------------------------------------------------------------------
  if args.collapse:
    print('Collapsing co-clustered mutations...', file=sys.stderr)
    # counter to create cluster ids from
    n_collapsed = 0
    # keep track of mutations in collapsed clusters
    collapsed_vars = {}
    print('1. Collapse mutations unique to each tree...', file=sys.stderr)
    # collapse unique vars in first tree
    true_vars_uniq = true_vars_set.difference(inf_vars_set)
    for tclust, tvars in true_clusters.items():
      v_uniq = true_vars_uniq & set(tvars)
      if len(v_uniq) > 1: # create new cluster for non-trivial set
        n_collapsed += 1
        id_clust = 'collapsed_{0}'.format(n_collapsed)
        collapsed_vars[id_clust] = list(v_uniq)
        true_snv_clust[(id_clust, len(v_uniq))] = tclust
        for v in v_uniq:
          true_snv_clust.pop(v)
    # collapse unique vars in second tree
    inf_vars_uniq = inf_vars_set.difference(true_vars_set)
    for iclust, ivars in inf_clusters.items():
      v_uniq = inf_vars_uniq & set(ivars)
      if len(v_uniq) > 1: # create new cluster for non-trivial set
        n_collapsed += 1
        id_clust = 'collapsed_{0}'.format(n_collapsed)
        collapsed_vars[id_clust] = list(v_uniq)
        inf_snv_clust[(id_clust, len(v_uniq))] = iclust
        for v in v_uniq:
          inf_snv_clust.pop(v)
    print('TRUE vars: {}'.format(len(true_vars_lst)), file=sys.stderr)
    print('INF vars:  {}'.format(len(inf_vars_lst)), file=sys.stderr)
    
    print('2. Collapse joint mutations...', file=sys.stderr)
    for tclust, tvars in true_clusters.items():
      for iclust, ivars in inf_clusters.items():
        isect = set(tvars) & set(ivars)
        if len(isect) > 1: # collapsing only makes sense for sets >= 2
          n_vars = len(isect)
          #id_var = '_'.join(list(isect))
          n_collapsed += 1
          id_clust = 'collapsed_{0}'.format(n_collapsed)
          collapsed_vars[id_clust] = list(isect)
          true_snv_clust[(id_clust, n_vars)] = tclust
          inf_snv_clust[(id_clust, n_vars)] = iclust
          for v in isect:
            true_snv_clust.pop(v)
            inf_snv_clust.pop(v)
    # update clusters using collapsed sets
    true_clusters.clear()
    for var, clust in true_snv_clust.items():
      if clust not in true_clusters:
        true_clusters[clust] = []
      true_clusters[clust].append(var)
    inf_clusters.clear()
    for var, clust in inf_snv_clust.items():
      if clust not in inf_clusters:
        inf_clusters[clust] = []
      inf_clusters[clust].append(var)
    true_vars_set = set(true_snv_clust.keys())
    inf_vars_set = set(inf_snv_clust.keys())
    print('TRUE vars: {}'.format(len(true_vars_set)), file=sys.stderr)
    print('INF vars:  {}'.format(len(inf_vars_set)), file=sys.stderr)
    
    # log collapsed clusters to file
    with open('collapsed_vars.csv', 'wt') as f:
      for id_clust, variants in collapsed_vars.items():
        for id_var, n in variants:
          f.write('{0},{1},{2}\n'.format(id_clust, id_var, n))

  # compile variant sets
  #-----------------------------------------------------------------------------
  snvs = true_vars_set
  snvs_isect = true_vars_set & inf_vars_set
  snvs_union = true_vars_set | inf_vars_set

  # calculate distance
  #-----------------------------------------------------------------------------
  if args.CASet_isect:
    tree_dist = CASet(snvs_isect, true_tree, inf_tree, true_clusters, inf_clusters, true_snv_clust, inf_snv_clust)
    print('CASet_isect: {:.4f}'.format(tree_dist))
  if args.CASet_union:
    tree_dist = CASet(snvs_union, true_tree, inf_tree, true_clusters, inf_clusters, true_snv_clust, inf_snv_clust)
    print('CASet_union: {:.4f}'.format(tree_dist))
  if args.DISC_isect:
    tree_dist = DISC(snvs_isect, true_tree, inf_tree, true_clusters, inf_clusters, true_snv_clust, inf_snv_clust)
    print('DISC_isect: {:.4f}'.format(tree_dist))
  if args.DISC_union:
    tree_dist = DISC(snvs_union, true_tree, inf_tree, true_clusters, inf_clusters, true_snv_clust, inf_snv_clust)
    print('DISC_union: {:.4f}'.format(tree_dist))


if __name__ == '__main__':
  args = parse_args()
  main(args)
