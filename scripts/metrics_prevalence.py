#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Calculate prevalence accuracy metrics.
#
# Input:
#   - TRUE clone prevalences (CSV: id_cluster,id_sample,freq)
#   - TRUE variant-to-cluster assignments (CSV: chrom_pos,id_cluster,...)
#   - INFERRED cluster prevalences (CSV: id_cluster,id_sample,freq)
#   - INFERRED variant-to-cluster assignment (CSV: chrom_pos,id_cluster,...)
# Output:
#   mean squared error between TRUE and INFERRED mutation CCFs
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-10-11
#------------------------------------------------------------------------------
import argparse
import numpy as np
import pandas as pd
from sklearn import metrics

def parse_args():
  parser = argparse.ArgumentParser(description='Calculate tree distance.')
  # positional arguments
  parser.add_argument('true_prev', type=argparse.FileType('r'), help='TRUE cluster prevalences (CSV: id_cluster,id_sample,freq).')
  parser.add_argument('true_snvs', type=argparse.FileType('r'), help='TRUE variant-to-cluster mapping (CSV: chrom_pos,id_cluster,...).')
  parser.add_argument('inf_prev', type=argparse.FileType('r'), help='INFERRED cluster prevalences (CSV: id_cluster,id_sample,freq).')
  parser.add_argument('inf_snvs', type=argparse.FileType('r'), help='INFERRED variant-to-cluster mapping (CSV: chrom_pos,id_cluster,...).')
  
  args = parser.parse_args()
  return args

def main(args):
  #print(vars(args))
  tru_snv = {}
  inf_snv = {}

  # read input assignments
  df_clust_tru = pd.read_csv(args.true_prev.name, usecols=[0,1,2], dtype={'id_cluster': str, 'id_sample': str, 'freq': np.float64})
  df_clust_inf = pd.read_csv(args.inf_prev.name, usecols=[0,1,2], dtype={'id_cluster': str, 'id_sample': str, 'freq': np.float64})
  df_snv_tru = pd.read_csv(args.true_snvs.name, usecols=[0,1], dtype={'chrom_pos': str, 'id_cluster': str})
  df_snv_inf = pd.read_csv(args.inf_snvs.name, usecols=[0,1], dtype={'chrom_pos': str, 'id_cluster': str})

  # assign cluster frequencies to mutations
  df_tru = pd.merge(df_snv_tru, df_clust_tru, on='id_cluster')
  df_tru = df_tru.rename(columns={'freq': 'freq_true'})
  df_inf = pd.merge(df_snv_inf, df_clust_inf, on='id_cluster')
  df_inf = df_inf.rename(columns={'freq': 'freq_inf'})
  
  # merge true and inferred data and calculate distance metric
  df = pd.merge(df_inf, df_tru, on=['chrom_pos', 'id_sample'])
  msq = metrics.mean_squared_error(df['freq_inf'], df['freq_true'])

  print('n_clust_true: {}'.format(len(df_clust_tru['id_cluster'].unique())))
  print('n_mut_true: {}'.format(len(df_snv_tru['chrom_pos'].unique())))
  print('n_clust_inf: {}'.format(len(df_clust_inf['id_cluster'].unique())))
  print('n_mut_inf: {}'.format(len(df_snv_inf['chrom_pos'].unique())))
  print('prev_msq: {:.4f}'.format(msq))

if __name__ == '__main__':
  args = parse_args()
  main(args)
