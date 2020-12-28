#!/usr/bin/env Rscript --vanilla
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Run SciClone for all input files in a given directory.
# Results are stored in same working directory.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-12-15
#------------------------------------------------------------------------------
library(sciClone)
library(reshape2)

# parse command line params
if(!exists("argv")) {
  argv <- commandArgs(trailingOnly = TRUE)
}
workdir <- as.character(argv[1])

# load result data
#-------------------------------------------------------------------------------
sc <- readRDS(file.path(workdir, "sc.rds"))

# extract variant-cluster assignments
df_mut <- sc@vafs.merged[,c("chr", "st", "cluster")]
df_mut$chrom_pos <- paste(df_mut$chr, df_mut$st, sep = "_")
df_mut$id_cluster <- df_mut$cluster
# extract cluster frequencies
m <- t(sc@clust[["cluster.means"]])
n_clust = nrow(m)
#rownames(m) = paste("cluster", 1:n_clust, sep="")
rownames(m) = 1:n_clust
colnames(m) = c(sc@sampleNames)
df_clust = melt(m, value.name = "freq", , varnames = c("id_cluster", "id_sample"))

# write output
#-------------------------------------------------------------------------------
#writeClusterTable(sc, file.path(workdir, 'result.clusters.tsv'))
#writeClusterSummaryTable(sc, file.path(workdir, 'result.clusters.summary'))
write.table(
  df_mut[,c('chrom_pos', 'id_cluster')],
  file.path(workdir, 'inf.snvs.csv'),
  sep = ',', 
  row.names = FALSE,
  quote = FALSE
)
write.table(
  df_clust,
  file.path(workdir, 'inf.clusters.csv'),
  sep = ',', 
  row.names = FALSE,
  quote = FALSE
)

