#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert an adjacency list plus a mutation list to a Newick mutation tree.
#
# Input:
#   - adjacency list (CSV: from,to)
#   - variant-to-cluster assignment (CSV: id_var,chrom,pos,id_clone)
# Output:
#   - Newick mutation tree
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-29
#------------------------------------------------------------------------------

import argparse
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Create Newick mutation tree.')
  parser.add_argument('edges', type=argparse.FileType('r'), help='Adjacency list (CSV: cluster_from,cluster_to).')
  parser.add_argument('clusters', type=argparse.FileType('r'), help='Variant-to-cluster mapping (CSV: id_var,id_cluster)')
  
  args = parser.parse_args()
  return args

def adjacency_list_to_tree(adj_list):
  '''Initialize ete3.Tree using a set of edges.'''
  nodes = {}
  for lbl_parent, lbl_child in adj_list:
    p = nodes[lbl_parent] if lbl_parent in nodes else Tree(name=lbl_parent)
    c = nodes[lbl_child] if lbl_child in nodes else Tree(name=lbl_child)
    p.add_child(c)
    nodes[lbl_parent] = p
    nodes[lbl_child] = c

  roots = [n for lbl, n in nodes.items() if n.is_root()]
  assert len(roots) == 1, "More than one root found."
  return roots[0]

def clusters_to_tree(edges, clusters):
  tree = adjacency_list_to_tree(edges)
  for node in tree.traverse('preorder'):
    variants = clusters[node.name] if node.name in clusters else []
    node.name = '{{{0}}}'.format(','.join(variants))
  
  return tree

def main(args):
  edges = []
  clusters = {}

  # read adjacency list
  for line in args.edges:
    row = line.strip().split(',')
    edges.append((row[0],row[1]))

  # read variant-to-cluster assignment
  for line in args.clusters:
    row = line.strip().split(',')
    id_var = row[0]
    id_clust = row[1]
    if id_clust in clusters:
      clusters[id_clust].append(id_var)
    else:
      clusters[id_clust] = [id_var]
  
  tree = clusters_to_tree(edges, clusters)
  # print Newick to stdout
  print(tree.write(format=1))

if __name__ == '__main__':
  args = parse_args()
  main(args)
