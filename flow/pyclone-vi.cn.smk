# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

OUTPFX = "pyclonevi_cn"

rule pyclonevi_cn_prep:
  input:
    "sequenza/snv.cn.csv",
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

rule pyclonevi_cn:
  input:
    "%s/input.snvs.tsv" % OUTPFX
  output:
    done = "%s/inf.done" % OUTPFX
  params:
    h5   = "%s/result.h5" % OUTPFX,
    tsv  = "%s/result.tsv" % OUTPFX
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
    pyclone-vi fit -i {input} -o {params.h5} -c 20 -d beta-binomial -r 10
    pyclone-vi write-results-file -i {params.h5} -o {params.tsv}

    touch {output.done}
    ) >{log} 2>&1
    """

rule pyclonevi_cn_post:
  input: "%s/inf.done" % OUTPFX
  output:
    prev="%s/inf.clusters.csv" % OUTPFX,
    snvs="%s/inf.snvs.csv" % OUTPFX
  params:
    in_h5  = "%s/result.h5" % OUTPFX,
    outdir = OUTPFX
  shell:
    """
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][pyclonevi_env]}
    set -eu

    python {config[scripts]}/pyclone-vi_write_result.py \
      --h5 {params.in_h5} --outdir {params.outdir}
    """

rule pyclonevi_cn_metrics:
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
