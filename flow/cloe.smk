# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

import os
from glob import glob
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.vcf.gz")
PWD = os.getcwd
OUTPFX = "cloe_cn2_seqz"

rule cloe_prep:
  input: "sequenza/snv.cn2.csv"
  output:
    alt = "%s/rc_alt.csv" % OUTPFX,
    tot = "%s/rc_tot.csv" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s_prep.log" % OUTPFX
  group: "cloe"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/snvs2cloe.py \
      --csv {input} \
      --normal "RN" \
      --outdir {params.dir}
    ) >{log} 2>&1
    """

rule cloe:
  input:
    alt="%s/rc_alt.csv" % OUTPFX,
    tot="%s/rc_tot.csv" % OUTPFX
  output:
    done = "%s/inf.done" % OUTPFX
  params:
    dir = OUTPFX,
    snv = "%s/inf.snvs_clones.csv" % OUTPFX,
    tre = "%s/inf.trees.csv" % OUTPFX,
    res = "%s/inf.clusters.csv" % OUTPFX
  log:
    "log/%s.log" % OUTPFX
  threads: 5
  group: "cloe"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][cloe_env]}
    set -eu

    Rscript {config[scripts]}/run_cloe.R {params.dir} {threads}
    touch {output.done}

    ) >{log} 2>&1
    """

rule cloe_post:
  input:
    done = "%s/inf.done" % OUTPFX
  output:
    snv = "%s/inf.snvs.csv" % OUTPFX
  params:
    snv = "%s/inf.snvs_clones.csv" % OUTPFX,
    tre = "%s/inf.trees.csv" % OUTPFX
  log:
    "log/%s_post.log" % OUTPFX
  shell:
    """
    time(
      python {config[scripts]}/tree_clones_to_clusters.py \
        {params.snv} {params.tre} > {output.snv}
    ) >{log} 2>&1
    """

rule cloe_metrics:
  input: "%s/inf.done" % OUTPFX
  output:
    clust = "%s/metrics_clustering.yml" % OUTPFX,
    prev  = "%s/metrics_prevalence.yml" % OUTPFX,
    phylo = "%s/metrics_phylogeny.yml" % OUTPFX
  params:
    inf_snv  = "%s/inf.snvs.csv" % OUTPFX,
    inf_prev = "%s/inf.clusters.csv" % OUTPFX,
    inf_tree = "%s/inf.trees.csv" % OUTPFX,
    tru_snv  = "sim/true.snvs.csv",
    tru_prev = "sim/true.clusters.csv",
    tru_tree = "sim/true.tree.csv"
  shell:
    """
    python {config[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {params.tru_snv} {params.inf_snv} > {output.clust}

    python {config[scripts]}/metrics_prevalence.py \
      {params.tru_prev} {params.tru_snv} {params.inf_prev} {params.inf_snv} \
    > {output.prev}
    
    python {config[scripts]}/metrics_phylogeny.py \
      --CASet_isect --CASet_union --DISC_isect --DISC_union \
      --collapse --ignore-homoplasy \
      {params.tru_tree} {params.tru_snv} {params.inf_tree} {params.inf_snv} \
    > {output.phylo}
    """
