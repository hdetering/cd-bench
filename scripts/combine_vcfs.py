#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Combine several single-sample VCFs into a multi-sample VCF.
#
# Requirements:
# - input VCFs have the same normal sample id
#
# Outputs:
# - VCF header of first infput VCF
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-04-21
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
from collections import OrderedDict
#import csv
import vcf, vcf.utils

def parse_args():
  parser = argparse.ArgumentParser(description='Create PyClone YAML input file from VCF and copy-number BED.')
  parser.add_argument('--vcf', required=True, action='append', 
    help='Input tumor-normal VCF file with somatic mutations and read counts.')
  parser.add_argument('--out', required=True, type=argparse.FileType('wt'), 
    help='Output multi-sample VCF file.')
  parser.add_argument('--bed', type=argparse.FileType(), 
    help='BED file with allele-specific copy number.')
  parser.add_argument('--info-to-format', action='append', 
    help='Map INFO fields to FORMAT fields (id_INFO:id_FORMAT,[...])')
  parser.add_argument('--nofilt', action='store_true', 
    help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  args_dict = vars(args)

  # two ways of inputting VCFs:
  # 1. filename
  # 2. sample_id:filename
  args_dict['samples'] = []
  args_dict['fn_vcf'] = []
  for vcf in args_dict['vcf']:
    if ':' in vcf:
      sid, fn = vcf.split(':')
      args_dict['samples'].append(sid)
    else:
      fn = vcf
    args_dict['fn_vcf'].append(fn)
  # make sure files exist
  for fn in args_dict['fn_vcf']:
    if not os.path.exists(fn):
      print('[ERROR] file does not exist: {}'.format(fn), file=sys.stderr)
      sys.exit(1)
  
  # parse INFO-to-FORMAT mapping
  args_dict['info-to-format'] = {}
  for fifo in args.info_to_format:
    fi, fo = fifo.split(':')
    args_dict['info-to-format'][fi] = fo

  return args_dict

def get_format_fields(rec, fields_info, fields_format=[], info_to_format={}, idx_sample=0):
  '''Extract relevant fields from a vcf.Record object.'''
  fields = {'INFO': {}, 'FORMAT': {}}
  # sanity checks
  assert type(rec) == vcf.Record

  for f in fields_format:
    if idx_sample < len(rec.samples) and f in rec.samples[idx_sample]:
      fields['FORMAT'][f] = rec.FORMAT[f]
  for f in fields_info:
    if f in rec.INFO:
      fields['INFO'][f] = rec.INFO[f]
  for fi, fo in info_to_format.items():
    if fi in rec.INFO:
      fields['FORMAT'][fo] = rec.INFO[fi]

  return fields

def get_out_header(rdr, info_to_format={}):
  hdr = rdr
  fmt = []
  for fi, fo in info_to_format.items():
    info = hdr.infos.pop(fi)
    hdr.formats[fo] = vcf.parser._Format(id=info.id, num=info.num, type=info.type, desc=info.desc)
    fmt.append(fo)
  if len(hdr._column_headers) == 8:
    hdr._column_headers.append('FORMAT')

  return hdr

def merge_field_allele_wise(merged_alleles, sample_alleles, sample_values):
  '''
  Update sample fields to apply to a merged set of alleles.
  
  Parameters:
    merged_alleles : joint list of alleles ([vcf.model._Substitution])
    sample_alleles : list of lists of alleles per sample ([[vcf.model._Substitution], ...])
    sample_values  : list of lists of value for each allele ()
  '''
  assert len(sample_alleles) == len(sample_values)
  new_vals = []
  for idx_sample in range(len(sample_alleles)):
    alleles = sample_alleles[idx_sample]
    values = sample_values[idx_sample]
    if alleles is None or values is None:
      new_vals.append(None)
      continue
    assert len(alleles) == len(values)
    v = [values[alleles.index(a)] if a in alleles else 0 for a in merged_alleles]
    new_vals.append(v)

  return new_vals

def merge_records(records, lbl_samples, infos, formats, info2format):
  '''Merge a list of vcf.Record objects into a single record.'''
  chrom = ""; pos = 0; id = ""; ref = ""; alt = []; info = {}; fmt = {}

  lst_chrom = []
  lst_pos = []
  lst_id = []
  lst_ref = []
  lst_alt = []
  #base_count = {lbl: {nuc: 0 for nuc in 'ACGT'} for lbl in lbl_samples}
  for i in range(len(records)):
    rec = records[i]
    if rec:
      lst_chrom.append(rec.CHROM)
      lst_pos.append(rec.POS)
      if rec.ID:
        lst_id.append(rec.ID)
      lst_ref.append(rec.REF)
      lst_alt += [str(sub) for sub in rec.ALT]

  # sanity checks
  assert len(set(lst_chrom)) == 1 # if this fails, coordinates don't match
  assert len(set(lst_pos)) == 1 # if this fails, coordinates don't match
  assert len(set(lst_ref)) == 1 # this can fail for non-SNVs, not clear how to proceed

  chrom = set(lst_chrom).pop()
  pos   = set(lst_pos).pop()
  id    = ','.join(list(set(lst_id)))
  ref   = set(lst_ref).pop()
  alt   = [vcf.model._Substitution(x) for x in set(lst_alt)]
  qual  = '.' # NOTE: try to calculate mean?
  filt  = 'PASS'

  merged_formats = OrderedDict()
  # parse info fields
  for fi, fo in info2format.items():
    field_vals = []
    if formats[fo].num == -1: # allele-wise field
      field_vals = merge_field_allele_wise(
        alt, [rec.ALT if rec else None for rec in records], [rec.INFO[fi] if rec else None for rec in records]
      )
    else: # sample-wise field
      field_vals = [rec.INFO[fi] if rec else None for rec in records]
    merged_formats[fo] = field_vals

  # contruct FORMAT field:
  fmt = ':'.join(merged_formats.keys())

  # construct merged Record
  idx_smp = {lbl_samples[i]: i for i in range(len(lbl_samples))}
  rec = vcf.model._Record(chrom, pos, id, ref, alt, qual, filt, info, fmt, idx_smp)
  CallDataType = vcf.model.make_calldata_tuple(merged_formats.keys())
  calldata_vals = list(zip(*merged_formats.values()))
  rec.samples = [
    vcf.model._Call(
      rec, 
      lbl_samples[i], 
      CallDataType(*calldata_vals[i]))
    for i in range(len(idx_smp))
  ]

  return rec

def main(args):
  vcf_in  = [vcf.Reader(filename=fn) for fn in args['fn_vcf']]
  # create header for output VCF
  hdr_out = get_out_header(vcf_in[0], args['info-to-format'])
  hdr_out.samples = args['samples']
  idx_smp = {args['samples'][i]: i for i in range(len(args['samples']))}
  hdr_out._sample_indexes = idx_smp
  vcf_out = vcf.Writer(args['out'], hdr_out)

  # if sample ids were provided by user, use those
  use_custom_sids = ( len(args['samples']) == len(args['fn_vcf']) )
  
  vars_comb = {} # {CHROM: {POS: {SAMPLE: var}}
  recs_out = []
  for recs in vcf.utils.walk_together(*(r for r in vcf_in)):
    r = merge_records(recs, args['samples'], hdr_out.infos, hdr_out.formats, args['info-to-format'])
    recs_out.append(r)
    #for i in range(len(recs)):
    #  sid = args['samples'][i] if use_custom_sids else ''
    #for var in rdr:
    #  if var.CHROM not in vars_comb:
    #    vars_comb[var.CHROM] = {}
    #  if var.POS not in vars_comb[var.CHROM]:
    #    vars_comb[var.CHROM][var.POS] = {}
    #  if use_custom_sids:
    #    # NOTE: maybe add sample index to "--vcf" param if multiple samples?
    #    vars_comb[var.CHROM][var.POS][sid] = ""
    vcf_out.write_record(r)

if __name__ == '__main__':
  args = parse_args()
  main(args)

