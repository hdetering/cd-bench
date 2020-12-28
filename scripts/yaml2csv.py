#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert YAML file to CSV format.
# 
# Expected YAML structure:
#   - { key1: value1_1, key2: value2_1, ...}
#   - { key1: value1_2, key2: value2_2, ...}
#
# Ouput CSV format:
#   key1,key2,...
#   value1_1,value2_1,...
#   value1_2,value2_2,...
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-22
#------------------------------------------------------------------------------

from __future__ import print_function
import argparse
import pandas as pd
import yaml
import yamlordereddictloader

def parse_args():
  parser = argparse.ArgumentParser(description='Convert YAML file to CSV format.')
  parser.add_argument('yaml', type=argparse.FileType(), help='YAML input file.')

  args = parser.parse_args()
  return args

def main(args):
  lst_yml = yaml.load(args.yaml, Loader=yamlordereddictloader.Loader)
  df = pd.DataFrame(lst_yml)
  #print(df.head())
  print(df.to_csv(index=False, columns=list(df.columns)))

if __name__ == '__main__':
  args = parse_args()
  main(args)
