#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert a Newick mutation tree to an adjacency list.
#
# Input:
#   - tree in Newick/NEXUS format
# Output:
#   - adjacency list in CSV format (from,to)
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-11
#------------------------------------------------------------------------------

import argparse
from nexus import NexusReader
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Create adjacency list from Newick tree.')
  grp_infile = parser.add_mutually_exclusive_group(required=True)
  grp_infile.add_argument('--newick', type=argparse.FileType('r'), help='File containing tree in Newick format.')
  grp_infile.add_argument('--nexus', type=argparse.FileType('r'), help='File containing tree in NEXUS format.')
  
  args = parser.parse_args()
  return args

def tree_to_adjacency_list(nwk, fmt=1):
  '''
  Extract edges from tree.
  
  Unnamed nodes will be labeled "node_x" (with x>=0).
  '''
  tree = Tree(nwk, format=fmt)
  pfx = 'node_' # prefix for labeling anonymous nodes
  i = 0 # index for labeling anonymous nodes
  edges = []
  for node in tree.traverse('preorder'):
    if not node.name:
      node.name = '{}{}'.format(pfx, i)
      i += 1
    for child in node.children:
      if not child.name:
        child.name = '{}{}'.format(pfx, i)
        i += 1
      edges.append((node.name, child.name))
  
  return edges

def main(args):
  newick = ''
  if args.newick:
    newick = args.newick.readline()
  elif args.nexus:
    fn_nexus = args.nexus.name
    nexus = NexusReader.from_file(fn_nexus)
    newick = nexus.trees.trees[0].split('=')[1].strip()

  edges = tree_to_adjacency_list(newick)
  # print adjacency list to stdout
  for v1, v2 in edges:
    print('{},{}'.format(v1, v2))

if __name__ == '__main__':
  args = parse_args()
  main(args)
