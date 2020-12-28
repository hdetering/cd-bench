#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Calculate clustering performance metrics.
#
# Input:
#   - TRUE variant-to-cluster assignments (CSV: chrom_pos,id_cluster,...)
#   - INFERRED variant-to-cluster assignment (CSV: chrom_pos,id_cluster,...)
# Output:
#   metrics value
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-21
#------------------------------------------------------------------------------
import argparse
from sklearn import metrics

def parse_args():
  parser = argparse.ArgumentParser(description='Calculate tree distance.')
  # positional arguments
  parser.add_argument('true_snvs', type=argparse.FileType('r'), help='TRUE variant-to-cluster mapping (CSV: chrom_pos,id_cluster,...).')
  parser.add_argument('inf_snvs', type=argparse.FileType('r'), help='INFERRED variant-to-cluster mapping (CSV: chrom_pos,id_cluster,...).')
  # optional arguments
  #grp_metric = parser.add_mutually_exclusive_group(required=True)
  parser.add_argument('--ARI', action='store_true', help='Calculate Adjusted Rand Index')
  parser.add_argument('--V-measure', action='store_true', help='Calculate V-measure (reports homogeneity and completeness, too)')
  
  args = parser.parse_args()
  return args

def main(args):
  #print(vars(args))
  tru_snv = {}
  inf_snv = {}

  # read input assignments
  hdr = args.true_snvs.readline()
  for line in args.true_snvs:
    cols = line.strip().split(',', 2)
    tru_snv[cols[0]] = cols[1]
  hdr = args.inf_snvs.readline()
  for line in args.inf_snvs:
    cols = line.strip().split(',', 2)
    inf_snv[cols[0]] = cols[1]

  # extract cluster assignments for vars in TRUE and INFERRED
  joint_snv = list(set(tru_snv.keys()) & set(inf_snv.keys()))
  tru_clust = [tru_snv[x] for x in joint_snv]
  inf_clust = [inf_snv[x] for x in joint_snv]

  print('n_clust_true: {}'.format(len(set(tru_snv.values()))))
  print('n_mut_true: {}'.format(len(tru_snv.keys())))
  print('n_clust_inf: {}'.format(len(set(inf_snv.values()))))
  print('n_mut_inf: {}'.format(len(inf_snv.keys())))

  if args.ARI:
    ari = metrics.adjusted_rand_score(tru_clust, inf_clust)
    print('ARI: {:.4f}'.format(ari))
  if args.V_measure:
    hom = metrics.homogeneity_score(tru_clust, inf_clust)
    com = metrics.completeness_score(tru_clust, inf_clust)
    v = metrics.v_measure_score(tru_clust, inf_clust)
    print('homogeneity: {:.4f}'.format(hom))
    print('completeness: {:.4f}'.format(com))
    print('V_measure: {:.4f}'.format(v))

if __name__ == '__main__':
  args = parse_args()
  main(args)
