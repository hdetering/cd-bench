# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

import os
from glob import glob
SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.vcf.gz")
OUTPFX = "phylowgs_cn_seqz"

rule phylowgs_prep:
  input:
   csv = "sequenza/snv.cn.csv",
  output:
    snv="%s/input.snvs.tsv" % OUTPFX,
    cnv="%s/input.cnvs.tsv" % OUTPFX
  params:
    dir=OUTPFX
  log: "log/%s_prep.log" % OUTPFX
  group: "phylowgs"
  shell:
    """
    time (
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/vcf2phylowgs.py \
      --csc {input.csv} \
    > {output.snv}
    # create dummy file for CNVs
    touch {output.cnv}
    #printf "cnv\ta\td\tssms\tphysical_cnvs\n" > {output.cnv} 
    ) >{log} 2>&1
    """

rule phylowgs:
  input:
    snv="%s/input.snvs.tsv" % OUTPFX,
    cnv="%s/input.cnvs.tsv" % OUTPFX
  output:
    done="inf.done",
  params:
    out_tre="%s/trees.zip" % OUTPFX,
    out_sum="%s/result.summ.json.gz" % OUTPFX,
    out_mut="%s/result.muts.json.gz" % OUTPFX,
    out_ass="%s/result.mutass.zip" % OUTPFX,
    outdir=OUTPFX
  log:
    "log/%s.log" % OUTPFX
  threads: 4
  group: "phylowgs"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][phylowgs_env]}
    set -eu

    python2 {config[tools][phylowgs_env]}/share/phylowgs/multievolve.py \
      --num-chains {threads} \
      --ssms {input.snv} \
      --cnvs {input.cnv} \
      --output-dir {params.outdir}

    python2 {config[tools][phylowgs_env]}/share/phylowgs/write_results.py \
      results {params.tre} {params.sum} {params.mut} {params.ass}
    ) >{log} 2>&1

    touch {output.done}
    """

rule phylowgs_post:
  input: "%s/inf.done"
  output:
    snvs="%s/inf.snvs.csv" % OUTPFX,
    prev="%s/inf.clusters.csv" % OUTPFX,
    tree="%s/inf.trees.csv" % OUTPFX
  params:
    in_sum="%s/result.summ.json.gz" % OUTPFX,
    in_mut="%s/result.muts.json.gz" % OUTPFX,
    in_ass="%s/result.mutass.zip" % OUTPFX
  shell:
    """
    python3 {config[scripts]}/phylowgs_parse_result \
      {params.in_mut} {params.in_sum} {params.in_ass}
    """

rule phylowgs_metrics:
  input:  
    inf_tree = "%s/inf.trees.csv" % OUTPFX,
    inf_snvs = "%s/inf.snvs.csv" % OUTPFX,
    inf_prev = "%s/inf.clusters.csv" % OUTPFX,
    true_tree = "sim/true.tree.csv",
    true_snvs = "sim/true.snvs.csv",
    true_prev = "sim/true.clusters.csv"
  output:
    clust = "%s/metrics_clustering.yml" % OUTPFX,
    tree  = "%s/metrics_phylogeny.yml" % OUTPFX,
    prev  = "%s/metrics_prevalence.yml" % OUTPFX
  shell:
    """
    set -x
    python3 {config[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {input.true_snvs} {input.inf_snvs} \
    | tee {output.clust}

    # calculate distance measures
    # NOTE: add '--collapse' if running too slowly
    python3 {config[scripts]}/metrics_phylogeny.py \
      --CASet_isect --CASet_union --DISC_isect --DISC_union --collapse\
      {input.true_tree} \
      {input.true_snvs} \
      {input.inf_tree} \
      {input.inf_snvs} \
    | tee {output.tree}

    python3 {config[scripts]}/metrics_prevalence.py \
      {input.true_prev} \
      {input.true_snvs} \
      {input.inf_prev} \
      {input.inf_snvs} \
    | tee {output.prev}
    """
