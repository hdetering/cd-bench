# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

import os
import bisect
from glob import glob
import yaml
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.vcf.gz")
PWD = os.getcwd
OUTPFX = "pyclonevi_cn2_seqz"

rule pyclonevi_prep:
  input:
    "sequenza/snv.cn2.csv",
  output:
    "%s/input.snvs.tsv" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s_prep.log" % OUTPFX
  group: "pyclone"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/snvs2pyclone-vi.py \
      --csv {input} \
      > {output}
    ) >{log} 2>&1
    """

rule pyclonevi:
  input:
    "%s/input.snvs.tsv" % OUTPFX
  output:
    h5   = "%s/result.h5" % OUTPFX,
    tsv  = "%s/result.tsv" % OUTPFX,
    done = "%s/inf.done" % OUTPFX 
  log:
    "log/%s.log" % OUTPFX
  threads: 1
  group: "pyclone"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][pyclonevi_env]}
    set -eu
    pyclone-vi fit -i {input} -o {output.h5} -c 20 -d beta-binomial -r 10
    pyclone-vi write-results-file -i {output.h5} -o {output.tsv}

    touch {output.done}
    ) >{log} 2>&1
    """

rule pyclonevi_post:
  input: "%s/result.h5" % OUTPFX,
  output:
    prev="%s/inf.clusters.csv" % OUTPFX,
    snvs="%s/inf.snvs.csv" % OUTPFX
  params:
    dir = OUTPFX
  shell:
    """
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][pyclonevi_env]}
    set -eu

    python3 {config[scripts]}/pyclone-vi_write_result.py \
      --h5 {input} --outdir {params.dir}
    """

rule pyclonevi_metrics:
  input:  
    inf_snvs = "%s/inf.snvs.csv" % OUTPFX,
    inf_prev = "%s/inf.clusters.csv" % OUTPFX,
    true_snvs = "sim/true.snvs.csv",
    true_prev = "sim/true.clusters.csv"
  output:
    clust = "%s/metrics_clustering.yml" % OUTPFX,
    prev  = "%s/metrics_prevalence.yml" % OUTPFX,
  shell:
    """
    set -x
    python3 {config[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {input.true_snvs} {input.inf_snvs} \
    | tee {output.clust}

    python3 {config[scripts]}/metrics_prevalence.py \
      {input.true_prev} \
      {input.true_snvs} \
      {input.inf_prev} \
      {input.inf_snvs} \
    | tee {output.prev}
    """
