#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert YAML file into CSV.
# Only simple "key: value" pairs are allowed. Keys become column headers.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-21
#------------------------------------------------------------------------------

# valid YAML format requires a colon (':')
/:/ {
  split($0, kv, ":");
  # make sure there's only one colon
  if (length(kv) != 2) {
    print("[ERROR] Invalid YAML format found:");
    print("[ERROR] "$0);
    print("[ERROR] Bailing out...");
    exit 1;
  }
  # extract key and value
  k = kv[1];
  v = kv[2];
  # strip off spaces
  gsub(/ /, "", k);
  gsub(/ /, "", v);
  # store for later
  keys[++i] = k;
  vals[i] = v;
}
# print key-value pairs transposed
END {
  # print header row
  for (j=1; j<i; j++) {
    printf(keys[j]",");
  }
  printf(keys[i]"\n");
  # print data row
  for (j=1; j<i; j++) {
    printf(vals[j]",");
  }
  printf(vals[i]"\n");
}
