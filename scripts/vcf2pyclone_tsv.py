#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert VCF to PyClone/MuClone YAML input.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-04-11
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import csv
import vcf
import yaml
#from sortedcontainers import SortedDict
from intervaltree import Interval, IntervalTree

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone YAML input file from VCF and copy-number BED.')
  parser.add_argument('--vcf', required=True, type=argparse.FileType(), help='VCF file with somatic mutations and read counts.')
  parser.add_argument('--bed', type=argparse.FileType(), help='BED file with allele-specific copy number.')
  parser.add_argument('--out', required=True, type=argparse.FileType('wt'), help='Output YAML file.')
  parser.add_argument('--sample', help='Which sample\'s variants to extract (default: detect from VCF header).')
  parser.add_argument('--use-info', action='store_true', help='Extract read counts from INFO field (overrides --sample)')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  return args


def get_sample_idx(lbl_sample, vcf_reader):
  '''Determine index of sample by label.'''
  idx_sample = 0
  if lbl_sample:
    if lbl_sample not in vcf_reader.samples:
      print('[WARN] sample {} not present in VCF file. Attempting to guess tumor sample from header.'.format(sample))
      lbl_sample = None
  # attempt to determine tumor sample from VCF header
  if not lbl_sample and 'tumor_sample' in vcf_reader.metadata:
    lbl_sample = vcf_reader.metadata['tumor_sample'][0]
  if lbl_sample in vcf_reader.samples:
    idx_sample = vcf_reader.samples.index(lbl_sample)
  
  return idx_sample

def parse_bed(fh_bed):
  '''Resd allele-specific copy number from BED file.

  Expected format: 
    CHROM START END A B
  Returns:
    {"chrom": IntervalTree((START:END): (A, B))} for fast binary lookups
  '''
  chr2cn_tree = {}
  csv_reader = csv.reader(fh_bed, delimiter='\t')
  for row in csv_reader:
    assert(len(row)>4)
    chrom = row[0]
    start = int(row[1])
    end   = int(row[2])
    cn_a  = round(float(row[3]))
    cn_b  = round(float(row[4]))
    if chrom not in chr2cn_tree:
      chr2cn_tree[chrom] = IntervalTree()
    chr2cn_tree[chrom][start:end] = (cn_a, cn_b)
  
  return chr2cn_tree

def get_copynum_at(chrom, pos, cn_index):
  '''Returns allele-specific copy number at given locus.'''
  #import pdb; pdb.set_trace()
  assert(chrom in cn_index)
  cn_a, cn_b = -1.0, -1.0
  set_itvl = cn_index[chrom][pos]
  if len(set_itvl) == 1:
    cn_a, cn_b = list(set_itvl)[0].data
  elif len(set_itvl) > 1:
    print('[ERROR] >1 intervals found for locus "%s:%d"' % (chrom, pos))

  return cn_a, cn_b
  
def main(fh_output, fh_vcf, fh_bed, use_info, sample=None, nofilt=False):
  '''Extract filter-passing biallelic SNVs from VCF file.'''
  muts = []

  cn_index = parse_bed(fh_bed)
  rdr = vcf.Reader(filename=fh_vcf.name)
  
  if use_info:
    print('[INFO] exporting read counts from "INFO" field.')
  else:
    idx_sample = get_sample_idx(sample, rdr)
    print('[INFO] exporting read counts for sample "{}".'.format(rdr.samples[idx_sample]))  
  #import pdb; pdb.set_trace()
  for rec in rdr:
    # check filtering criteria
    if not nofilt and rec.FILTER and len(rec.FILTER) > 0:
      continue
    for idx_alt in range(len(rec.ALT)):
      # make sure variant is SNV
      if rec.ALT[idx_alt].type != 'SNV':
        continue
      # use INFO field or sample FORMAT fields
      if use_info:
        rc_ref = rec.INFO['DP']
        rc_alt = rec.INFO['AC'][idx_alt]
      else:
        call = rec.samples[idx_sample]
        rc_ref = call.data.AD[0]
        rc_alt = call.data.AD[idx_alt]
      # determine copy number at mutation locus
      cn_a, cn_b = get_copynum_at(rec.CHROM, rec.POS, cn_index)
      # at least one copy must be present (otherwise PyClone will complain)
      if cn_a + cn_b == 0:
        continue
      muts.append([
        "%s:%d_%d" % (rec.CHROM, rec.POS, idx_alt),
        rc_ref, rc_alt, 2, min([cn_a, cn_b]), max([cn_a, cn_b])
      ])
  header = ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']
  fh_output.write('\t'.join(header) + '\n')
  for row in muts:
    fh_output.write('\t'.join([str(x) for x in row]) + '\n')

if __name__ == '__main__':
  args = parse_args()
  main(args.out, args.vcf, args.bed, args.use_info, args.sample, args.nofilt)
