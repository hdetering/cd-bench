#!/usr/bin/env Rscript --vanilla
# vim: syntax=r tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Run Cloe for input files in a given directory.
# Results are stored in same working directory.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-31
#------------------------------------------------------------------------------
library(parallel)
library(cloe)

# parse command line params
if(!exists("argv")) {
  argv <- commandArgs(trailingOnly = TRUE)
}
workdir <- as.character(argv[1])
# number of cores to use
num_threads <- ifelse(length(argv)>1, argv[2], 1)
print(paste("Running in parallel using", num_threads, "cores."))

# run Cloe for this range of values of K
ks = 2:10
# output file names
fn_snvs <- file.path(workdir, 'inf.snvs_clones.csv')
fn_trees <- file.path(workdir, 'inf.trees.csv')
fn_clusters <- file.path(workdir, 'inf.clusters.csv')

# load input data
#-------------------------------------------------------------------------------
fn_reads <- file.path(workdir, 'rc_alt.csv')
fn_depths <- file.path(workdir, 'rc_tot.csv')
reads <- as.matrix(read.table(fn_reads, sep=',', header=T, row.names=1))
depths <- as.matrix(read.table(fn_depths, sep=',', header=T, row.names=1))

# create input object
#-------------------------------------------------------------------------------
ci <- cloe_input$new(reads, depths)
# cluster mutations
cic <- crp(ci)
saveRDS(cic, file.path(workdir, 'cic.rds'))

# run Cloe
#-------------------------------------------------------------------------------
# run MCMC chains
#cms <- lapply(ks, function(k) {
cms <- mclapply(ks, mc.cores=num_threads, FUN=function(k) {
  sampler(input=cic, iterations=1e5, K=k)
})
# write MCMC result to file (backup in case later steps fail)
saveRDS(cms, file.path(workdir, 'cms.rds'))

# get best sets of parameters
css <- mclapply(cms, mc.cores=num_threads, FUN=function(cmx) {
  summarise(cm=cmx, burn=0.5, thin=2, solutions=5L)
})

# select model
#-------------------------------------------------------------------------------
top_cs <- select_model(l=css, solutions=1L, plot=T)

# write results to output files
#-------------------------------------------------------------------------------
library(tidyverse)
# trees
if (!top_cs[[1]]$get_phylogenies()[[1]] %>% is.matrix()) {
top_cs[[1]]$get_phylogenies()[[1]] %>%
  matrix(nrow=1) %>%
  as_tibble() %>%
  rename(from = V1, to = V2) %>%
  write_csv(fn_trees)
} else {
top_cs[[1]]$get_phylogenies()[[1]] %>%
  as_tibble() %>% 
  rename(from = V1, to = V2) %>% 
  write_csv(fn_trees)
}
# clusters
df_clust <- top_cs[[1]]$get_clonal_fractions()[[1]] %>%
  as_tibble(rownames = 'id_cluster')
df_clust %>% 
  gather(-id_cluster, key = 'id_sample', value = 'freq') %>%
  write_csv(fn_clusters)
# variants
mm <- map_back(top_cs[[1]]$get_genotypes()[[1]], cic$get_clusters())
rownames(mm) <- ci$get_names()
colnames(mm) <- df_clust$id_cluster
mm %>% as_tibble(rownames = 'chrom_pos') %>%
  gather(-chrom_pos, key = 'id_cluster', value = 'present') %>%
  dplyr::filter(present == 1) %>%
  select(chrom_pos, id_cluster) %>%
  write_csv(fn_snvs)

