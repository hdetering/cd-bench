# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

import os
from glob import glob
SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.vcf.gz")
PWD = os.getcwd
OUTPFX = "sciclone_cn2_seqz"

rule sciclone_prep:
  input:
    csv = "sequenza/snv.cn2.csv",
  output:
    ["%s/%s.vaf.csv" % (OUTPFX, x) for x in SAMPLES]
  params:
    dir = OUTPFX
  log: "log/%s_prep.log" % OUTPFX
  group: "sciclone"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/snvs2sciclone.py \
      --csv {input.csv} \
      --outdir {params.dir}
      #--normal "RN"
    ) >{log} 2>&1
    """

rule sciclone:
  input:
    ["%s/%s.vaf.csv" % (OUTPFX, x) for x in SAMPLES]
  output:
    done = "%s/inf.done" % OUTPFX
  params:
    prev = "%s/inf.clusters.tsv" % OUTPFX,
    snvs = "%s/inf.snvs.csv" % OUTPFX,
    workdir = OUTPFX
  log:
    "log/%s.log" % OUTPFX
  threads: 1
  group: "sciclone"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][sciclone_env]}
    set -eu

    Rscript {config[scripts]}/run_sciclone.R {params.workdir}
    touch {output.done}

    ) >{log} 2>&1
    """

rule sciclone_metrics:
  input:
    "%s/inf.done" % OUTPFX
  params:
    inf_snvs = "%s/inf.snvs.csv" % OUTPFX,
    inf_prev = "%s/inf.clusters.csv" % OUTPFX,
    tru_snvs = "sim/true.snvs.csv",
    tru_prev = "sim/true.clusters.csv"
  output:
    clust = "%s/metrics_clustering.yml" % OUTPFX,
    prev  = "%s/metrics_prevalence.yml" % OUTPFX
  shell:
    """
    python {config[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {params.tru_snvs} {params.inf_snvs} \
    | tee {output.clust}

    python3 {config[scripts]}/metrics_prevalence.py \
      {params.tru_prev} \
      {params.tru_snvs} \
      {params.inf_prev} \
      {params.inf_snvs} \
    | tee {output.prev}
    """
