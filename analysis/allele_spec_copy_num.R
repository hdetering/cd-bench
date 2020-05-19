library(tidyverse)
library(tidygenomics)

data_dir <- '/home/harry/code/cd-bench/data'

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------

# find replicates
df_rep <- tibble( path = list.dirs(data_dir, full.names = TRUE, recursive = FALSE) ) %>%
  mutate( id_rep = str_extract(basename(path), '[^-]+') )

# read BED files
fn_bed <- list.files( 
  path = data_dir, 
  pattern = 'R\\d+.cn.bed', 
  recursive = TRUE,
  full.names = TRUE )
col_names <- c('chr', 'start', 'end', 'A', 'B')

df_bed <- tibble( fn = fn_bed ) %>% head() %>%
  extract( fn, c('id_rep'), '([^/-]+)-[^/]+/sim/bed', remove = FALSE ) %>%
  extract( fn, c('sample'), '(R[[:alnum:]]+).cn.bed$', remove = FALSE ) %>%
  mutate( bed = map(fn, ~ read_tsv(., col_names = col_names, col_types = 'ciidd')) ) %>%
  select( id_rep, sample, bed ) %>%
  unnest( cols = c(bed) )

df_seg <- df_rep %>% head(2) %>%
  mutate( seg = map(path, ~ { 
    read_tsv(file.path(., 'sim', 'segments.tsv'), 
             col_names = c('chr', 'start', 'end', 'id_clone', 'allele', 'id_seg'),
             col_types = 'ciiccc')
    }) ) %>%
  select( -path ) %>%
  unnest( cols = c(seg) )

df_var <- df_rep %>% head(2) %>%
  mutate( var = map(path, ~ { 
    read_tsv(file.path(., 'sim', 'segment_vars.tsv'), 
             col_names = c('id_seg', 'id_var', 'chrom', 'pos', 'ref', 'alt'),
             col_types = 'cccicc')
  }) ) %>%
  select( -path ) %>%
  unnest( cols = c(var) )

df_seg_var <- df_seg %>%
  inner_join( df_var, by = c('id_rep', 'id_seg') )

# sanity check: 
# each somatic variant should only exist on one allele
-------------------------------------------------------
# should be empty
df_seg_var %>% dplyr::filter( str_detect(id_var, '^s') ) %>%
  distinct( id_rep, id_var, allele ) %>%
  group_by( id_rep, id_var ) %>%
  tally() %>% 
  dplyr::filter( n > 1 )

#-------------------------------------------------------------------------------
# ALLELE-SPECIFIC COPY NUMBER
#-------------------------------------------------------------------------------

df_var_clone_cn <- df_seg_var %>%
  group_by( id_rep, id_var, ref, alt, allele ) %>%
  tally()

# get total copy number at variant loci
#---------------------------------------
df_var_nest <- df_var %>% 
  select( id_rep, id_var, chr = chrom, start = pos, end = pos ) %>%
  nest( variants = c(id_var, chr, start, end) )
df_seg_nest <- df_seg %>% 
  select( id_rep, id_clone, chr, start, end, allele ) %>%
  nest( segments = c(id_clone, chr, start, end, allele) )
df_seg_var_nest <- df_seg_nest %>%
  inner_join( df_var_nest, by = c('id_rep') )

df <- pmap_dfr( df_seg_var_nest, function(id_rep, segments, variants) {
  genome_intersect( variants, segments, by = c('chr', 'start', 'end'), mode = 'both' ) %>%
    mutate( id_rep = id_rep )
} )
df %>% group_by( id_rep, id_var, id_clone, allele ) %>%
  tally() %>% head()
