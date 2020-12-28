#!/usr/bin/env Rscript --vanilla
# vim: syntax=r tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Run Canopy for input files in a given directory.
# Results are stored in same working directory.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-09-13
#------------------------------------------------------------------------------
library(Canopy)
library(future)

# parse command line params
if(!exists("argv")) {
  argv <- commandArgs(trailingOnly = TRUE)
}
workdir <- as.character(argv[1])
# number of cores to use
num_threads <- ifelse(length(argv)>1, argv[2], 1)
print(paste("Running in parallel using", num_threads, "cores."))

# run Canopy for this range of values of K
ks = 2:10
prj = 'canopy'  # project name
num_chains = 10 # number of MCMC chains to run
fn_postmcmc = file.path(workdir, paste0(prj, '_postmcmc_image.rda')) # file with MCMC results
# output file names
fn_snvs <- file.path(workdir, 'canopy.snvs.csv')
fn_trees <- file.path(workdir, 'canopy.trees.csv')
fn_clusters <- file.path(workdir, 'canopy.clusters.csv')

# only run Canopy if it has not run before
if (!file.exists(fn_postmcmc)) {
# load input data
#-------------------------------------------------------------------------------
fn_reads <- file.path(workdir, 'rc_alt.csv')
fn_depths <- file.path(workdir, 'rc_tot.csv')
reads <- as.matrix(read.table(fn_reads, sep=',', header=T, row.names=1))
depths <- as.matrix(read.table(fn_depths, sep=',', header=T, row.names=1))
# create dummy data structures for CNAs
WM = Wm = matrix(
  rep(1, dim(reads)[2]),
  nrow = 1, 
  dimnames = list(c('cna1'), colnames(reads)))
Y = matrix(
  c(rep(c(1,0), dim(reads)[1])),
  ncol = 2,
  byrow = T,
  dimnames = list(rownames(reads), c('non-cna', 'cna1')))
eM = em = 0.02

# run Canopy
#-------------------------------------------------------------------------------
sampchain <- canopy.sample.parallel(
  R = reads, 
  X = depths,
  WM = WM, Wm = Wm,
  epsilonM = eM, epsilonm = em,
  C = NULL, 
  Y = Y,
  K = ks,
  numchain = num_chains,
  max.simrun = 10000000,
  min.simrun =  1000000,
  writeskip = 200,
  projectname = prj,
  cell.line = FALSE,
  plot.likelihood = TRUE
)
save.image(file = fn_postmcmc, compress = 'xz')
} else {
  print("Canopy was run before. Loading results from previous run...")
}

# select model (using BIC) to determine number of subclones
#-------------------------------------------------------------------------------
load(fn_postmcmc) # load results from MCMC runs
burnin = 10 # NOTE: burnin applied after 'writeskip'
thin = 5    # NOTE: thinning applied in addition to 'writeskip'
bic = canopy.BIC(
  sampchain = sampchain,
  projectname = prj,
  K = ks,
  numchain = num_chains,
  burnin = burnin,
  thin = thin,
  pdf = TRUE
)
optK = ks[which.max(bic)]

# posterior tree evaluation
#-------------------------------------------------------------------------------
post = canopy.post(
  sampchain = sampchain,
  projectname = prj,
  K = ks,
  numchain = num_chains,
  burnin = burnin,
  thin = thin,
  optK = optK,
  post.config.cutoff = 0.05
)
samptreethin = post[[1]]     # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]] # likelihoods of trees in samptree
config = post[[3]]
config.summary = post[[4]]
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space
# note: if modes of posterior probabilities aren't obvious, run sampling longer.

# tree output and plot
#-------------------------------------------------------------------------------
# choose the configuration with the highest posterior likelihood
config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree = canopy.output(post, config.i, C = NULL)
pdf.name = file.path(workdir, paste0(prj, '_config_highest_likelihood.pdf'))
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, pdf = FALSE)

# write results to output files
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggtree)

# extract tree edges, assign labels for internal nodes and tip clones
# NOTE: 'clone_1' is healthy cell clone (interpret CCF as contamination)
tree_data <- fortify(output.tree) %>%
  mutate(idx = cumsum(is.na(label))) %>%
  mutate(id_node = ifelse(is.na(label), paste0('node_', idx), paste0('clone_', label)))

# trees
e <- output.tree$edge
colnames(e) <- c('idx_from', 'idx_to')
e %>% as_tibble() %>%
  inner_join(tree_data %>% select(idx_from = node, id_node)) %>% rename(from = id_node) %>%
  inner_join(tree_data %>% select(idx_to = node, id_node)) %>% rename(to = id_node) %>%
  select(from, to) %>%
  write_csv(fn_trees)
# clusters
output.tree$P %>% as_tibble(rownames = 'id_cluster') %>%
  gather(-id_cluster, key = 'id_sample', value = 'freq') %>%
  write_csv(fn_clusters)
# variants
output.tree$clonalmut %>%
  map_df(enframe, .id = 'idx_clone') %>% unnest() %>%
  dplyr::filter(value != 'None') %>%
  mutate(id_cluster = paste0('clone_', idx_clone)) %>%
  select(chrom_pos = value, id_cluster) %>%
  write_csv(fn_snvs)
