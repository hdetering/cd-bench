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
import vcf
import yaml

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
  

def main(fh_output, fh_vcf, fh_bed, use_info, sample=None, nofilt=False):
  '''Extract filter-passing biallelic SNVs from VCF file.'''
  muts = []
  rdr = vcf.Reader(filename=fh_vcf.name)
  
  if use_info:
    print('[INFO] exporting read counts from "INFO" field.')
  else:
    idx_sample = get_sample_idx(sample, rdr)
    print('[INFO] exporting read counts for sample "{}".'.format(rdr.samples[idx_sample]))  
  import pdb; pdb.set_trace()
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
      muts.append({
        'id': '%s:%d_%d' % (rec.CHROM, rec.POS, idx_alt),
        'ref_counts': rc_ref,
        'var_counts': rc_alt,
        #'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight':1} for x in ['AB', 'BB']]
        'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight':1} for x in ['AB']]
      })
  data_out = {'mutations': muts}
  fh_output.write(yaml.dump(data_out))

if __name__ == '__main__':
  args = parse_args()
  main(args.out, args.vcf, args.bed, args.use_info, args.sample, args.nofilt)
