#!/usr/bin/env Rscript --vanilla
# vim: syntax=r tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Run Sequenza for given input file.
# Results are stored in same working directory.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-10-27
#------------------------------------------------------------------------------
library(sequenza)

# parse command line params
if(!exists("argv")) {
  argv <- commandArgs(trailingOnly = TRUE)
}
infile <- as.character(argv[1])
infile <- normalizePath(infile)

# extract sample id from filenames
id_sample <- gsub('(.+)\\.seqz\\.gz', '\\1', basename(infile))

seqz <- sequenza.extract(infile, window = 1000, verbose = FALSE)
cp <- sequenza.fit(seqz, segment.filter=1000)
sequenza.results(sequenza.extract = seqz,
    cp.table = cp, sample.id = id_sample,
    out.dir = dirname(infile))
