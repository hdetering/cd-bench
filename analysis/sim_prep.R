require(phytools)
require(ape)
require(ggtree)
require(tidyverse)
require(yaml)

sim_dir = "/home/harry/code/cd-bench/data/100"
sim_dir = "/home/harry/Dropbox/PhD/cd-bench/sims/100"
ana_dir = "/home/harry/code/cd-bench/analysis/100"

df_reps <- tibble(
  id       = character(),
  nclones  = numeric(),
  nsamples = numeric(),
  cvg      = numeric(),
  fn_tree  = character()
)
df_prev <- tibble(
  id        = character(),
  id_sample = character(),
  id_clone  = character(),
  prev      = numeric()
)
df_trees <- tibble(
  id            = character(),
  node          = integer(),
  parent        = integer(),
  branch.length = numeric(),
  x = numeric(), y = numeric(),
  label         = character(),
  isTip         = logical(),
  branch        = numeric(),
  angle         = numeric()
)

repdirs <- dir(sim_dir, include.dirs=F)
for (rep_dir in repdirs) {
  # read meta info about replicate
  fn_meta <- file.path(sim_dir, rep_dir, 'meta.yml')
  meta <- read_yaml(fn_meta)
  # read config file for replicate
  fn_conf <- file.path(sim_dir, rep_dir, 'config.yml')
  conf <- read_yaml(fn_conf)
  # construct clone tree file name
  fn_tree <- file.path(sim_dir, rep_dir, 'clone_tree.nwk')
  tree <- phytools::read.newick(fn_tree)

  # marshal relevant info for replicate
  df_reps <- df_reps %>% rbind(tibble(
      id       = rep_dir,
      nclones  = meta$nclones,
      nsamples = meta$nsamples,
      cvg      = conf$`seq-coverage`,
      fn_tree  = fn_tree))

  fn_prev <- file.path(sim_dir, rep_dir, 'clone_prev.csv')
  prev <- read.table(fn_prev, sep = ',', header = T, row.names = 1)
  prev <- prev %>% rownames_to_column(var = 'id_sample') %>% mutate(id = rep_dir)
  prev <- prev %>% gather(id_clone, prev, -id, -id_sample)
  df_prev <- df_prev %>% rbind(prev)

  treedata <- fortify(tree)
  treedata$id <- rep_dir
  df_trees <- df_trees %>% rbind(treedata)
}

# create index column for replicates within groups
df_reps <- df_reps %>% 
  dplyr::group_by(cvg) %>% 
  mutate(idx_rep = row_number()) %>% 
  ungroup()
df_reps %>% saveRDS(file.path(ana_dir, 'df_reps.rds'))
df_prev %>% saveRDS(file.path(ana_dir, 'df_prev.rds'))
df_trees %>% saveRDS(file.path(ana_dir, 'df_trees.rds'))

df <- df_prev %>% inner_join(df_reps %>% select(id, idx_rep, ttype, cvg), by='id')
p <- ggplot(df, aes(x = id_clone, y = id_sample)) +
  geom_tile(aes(fill = prev)) +
  scale_fill_gradient(low = 'white', high = 'red', name = 'prevalence') +
  facet_grid(cvg+idx_rep~ttype) + theme_bw()
ggsave("sim_prep.prev.pdf", plot = p, device = 'pdf', width = 8, height = 40)

df <- df_trees %>% inner_join(df_reps %>% select(id, idx_rep, ttype, cvg), by='id')
x_max <- max(df$x)
y_max <- max(df$y)
df <- df %>% group_by(id) %>%
  mutate(x=x+.5*(x_max-max(x)), y=y+.5*(y_max-max(y))) %>%
  ungroup()
p <- ggtree(df) + geom_label(aes(label=label), size=2) +
  xlim(0, (x_max+.1)) + ylim(.9, (y_max+.1)) +
  facet_grid(cvg+idx_rep~ttype) + theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("sim_prep.tree.pdf", plot = p, device = 'pdf', width = 8, height = 40)

# VAF spectra
#-------------------------------------------------------------------------------
df_reps <- readRDS(file.path(ana_dir, 'df_reps.rds'))
df_prev <- readRDS(file.path(ana_dir, 'df_prev.rds'))
# load read count data
fn_var <- file.path(ana_dir, 'snv.cn2.strict.csv')
df_var <- read_csv(fn_var, col_types = 'ccicccii') %>%
  mutate(vaf = ifelse(rc_tot > 0, rc_alt/rc_tot, 0.0))

p <- df_var %>% dplyr::filter(vaf > 0) %>%
  ggplot(aes(x = vaf)) +
  geom_histogram(bins = 50) +
  geom_vline(data = df_prev %>% mutate(id_rep = id), aes(xintercept = 0.5*prev, color = id_clone)) +
  facet_grid(id_rep ~ id_sample, scales = 'free_y') +
  theme(strip.text.y = element_text(angle = 0))
ggsave('sim.vaf_spectra.pdf', plot = p, device = 'pdf', width = 24, height = 100, limitsize = FALSE)
