# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# author : Harald Detering (harald.detering@gmail.com)
# date   : 21.01.2019
#------------------------------------------------------------------------------

from __future__ import division, print_function
import os, sys, argparse
import re

def parse_args():
  parser = argparse.ArgumentParser(description='Convert CNT-MD output file to BED-formatted file.')
  parser.add_argument('cntmd_input', 
                      metavar='CNTMD-input', 
                      type=argparse.FileType('r'), 
                      help='Input file used with CNT-MD (contains bp coordinates).')
  parser.add_argument('cntmd_output', 
                      metavar='CNTMD-output', 
                      type=argparse.FileType('r'), 
                      help='Output text file from CNT-MD (contains CN states).')
  args = parser.parse_args()

  return args

def parse_cntmd_input(file_handle):
  # extract CN profiles, CN profile tree and profile mixture from
  regex = r'''#PARAMS
(?P<nchrom>\d+)( #.*)?
(?P<nsamples>\d+) (#.*)?
(?P<nsegs>(?:\d+ )+)(#.*)?
#SAMPLES (?P<segments>[^\n]+)
(?P<samples>.+)'''
  m = re.match(regex, file_handle.read(), flags=re.MULTILINE | re.DOTALL)
  if not m:
    print('[ERROR] CNT-MD input file does not match expected format.', file=sys.stderr)
    sys.exit(1)
  
  nchrom = int(m.group('nchrom'))
  nsamples = int(m.group('nsamples'))
  nsegs = [int(x) for x in m.group('nsegs').strip().split(' ')]
  segments = []
  for x in  m.group('segments').strip().split(' | '):
    chrom, str_segs = x.strip().split(' : ')
    segments.append([chrom, [x.split(',') for x in str_segs.strip().split(' ')]])
  print("chromosomes: %s" % segments)

def parse_cntmd_output(file_handle):
  # extract CN profiles, CN profile tree and profile mixture from
  regex = r'''#PARAMS
(?P<nchrom>\d+)( #.*)?
(?P<nleaves>\d+)( #.*)?
(?P<nsegs>(?:\d+ )+)(#.*)?
#PROFILES
(?P<profiles>.+)
#EDGES
(?P<edges>.+)
#EVENTS
(?P<events>.+)
#MIX\-PARAMS
(?P<mix_params>.+)
#PROPORTIONS
(?P<proportions>.+)
#SAMPLES
(?P<samples>.+)'''
  m = re.match(regex, file_handle.read(), flags=re.MULTILINE | re.DOTALL)
  if not m:
    print('[ERROR] CNT-MD output file does not match expected format.', file=sys.stderr)
    sys.exit(1)
  
  nchrom = int(m.group('nchrom'))
  nleaves = int(m.group('nleaves'))
  nsegs = [int(x) for x in m.group('nsegs').strip().split(' ')]
  # PROFILES
  profiles = {}
  for x in m.group('profiles').strip().split('\n'):
    y, z = x.strip().split(' : ')
    profiles[int(y)] = z
  nprofiles = len(profiles.keys())
  # EDGES
  v1 = set()
  v2 = set()
  for x in m.group('edges').strip().split('\n'):
    a, b = x.strip().split(' -> ')
    v1.add(a)
    v2.add(b)
  vertices = v1 | v2
  leaves = vertices - v1
  assert( len(vertices) == nprofiles )
  assert( len(leaves) == nleaves )

  proportions = []
  for x in m.group('proportions').strip().split('\n'):
    proportions.append(x.strip().split(' '))

def main(fh_input, fh_output):
  print("Parsing Input file...")
  parse_cntmd_input(fh_input)
  print("Parsing Output file...")
  parse_cntmd_output(fh_output)

if __name__ == '__main__':
  args = parse_args()
  main(args.cntmd_input, args.cntmd_output)
