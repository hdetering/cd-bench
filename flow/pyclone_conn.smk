# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

import os
import bisect
from glob import glob
import yaml
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.vcf.gz")
PWD = os.getcwd
OUTDIR = "pyclone_conn_cn2_seqz"

rule pyclone_conn_prep:
  input:
    csv="sequenza/snv.cn2.csv",
  output:
    ["%s/%s.tsv" % (OUTDIR, x) for x in SAMPLES]
  params:
    dir=OUTDIR
  log: "log/%s_prep.log" % OUTDIR
  group: "pyclone"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/snvs2pyclone.py \
      --csv {input.csv} \
      --outdir {params.dir} \
      #--normal "RN"
    ) >{log} 2>&1
    """

rule pyclone_conn_yaml:
  input:  "%s/{sample}.tsv" % OUTDIR
  output: "%s/{sample}.yaml" % OUTDIR
  log:    "log/%s_yaml.{sample}.log" % OUTDIR
  group: "pyclone"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate pyclone
    set -eu
    PyClone build_mutations_file \
      --in_file {input} \
      --out_file {output} \
      --prior total_copy_number
    ) >{log} 2>&1
    """

rule pyclone_conn_conf:
  input:
    yml=["%s/%s.yaml" % (OUTDIR, x) for x in SAMPLES]
  output:
    cfg="%s/config.yaml" % OUTDIR
  group: "pyclone"
  run:
    import os, yaml
    cfg = {
      # Specifies working directory for analysis. All paths in the rest of the file are relative to this.
      'working_dir': os.path.join(os.getcwd(), OUTDIR),
      # Where the trace (output) from the PyClone MCMC analysis will be written.
      'trace_dir'  : 'trace',
      # Specifies which density will be used to model read counts. Most people will want pyclone_beta_binomial or pyclone_binomial
      #'density'    : 'pyclone_beta_binomial',
      'density'    : 'pyclone_binomial',
      # Number of iterations of the MCMC chain.
      #'num_iters'  : 10000,
      'num_iters'  : 25000,
      # How to initialize clusters (connected: all muts in single cluster, disconnected: each mut in own cluster)
      'init_method' : 'connected',

      # Specifies parameters in Beta base measure for DP. Most people will want the values below.
      'base_measure_params': {
        'alpha' : 1,
        'beta'  : 1
      },

      # Specifies initial values and prior parameters for the prior on the concentration (alpha) parameter in the DP. If the prior node is not set the concentration will not be estimated and the specified value will be used.
      'concentration': {
        # Initial value if prior is set, or fixed value otherwise for concentration parameter.
        'value' : 1.0,
        # Specifies the parameters in the Gamma prior over the concentration parameter.
        'prior' : {
          'shape' : 1.0,
          'rate'  : 0.001
        }
      },
      
       # Beta-Binomial precision (alpha + beta) prior
      'beta_binomial_precision_params': {
        # Starting value
        'value': 1000,
        # Parameters for Gamma prior distribution
        'prior': {
          'shape': 1.0,
          'rate': 0.0001
        },
        # Precision of Gamma proposal function for MH step
        'proposal': {
          'precision': 0.01
        }
      },

      # Contains one or more sub-entries which specify details about the samples used in the analysis.
      'samples': { x: {
        'mutations_file': '%s.yaml' % x,
        'tumour_content': { 'value': 1.0 },
        'error_rate': 0.01
        } for x in SAMPLES}
    }
    with open(output.cfg, 'wt') as f_out:
      f_out.write( yaml.dump(cfg) )

rule pyclone_conn:
  input:
    cfg="%s/config.yaml" % OUTDIR
  output:
    loci="%s/result.loci.tsv" % OUTDIR,
    clust="%s/result.clusters.tsv" % OUTDIR
  params:
    done = "%s/inf.done" % OUTDIR
  log:
    "log/%s.log" % OUTDIR
  threads: 1
  group: "pyclone"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][pyclone_env]}
    set -eu
    PyClone run_analysis --config_file {input.cfg} 
    PyClone build_table --config_file {input.cfg} --out_file {output.loci} --table_type loci
    PyClone build_table --config_file {input.cfg} --out_file {output.clust} --table_type cluster

    touch {params.done}
    ) >{log} 2>&1
    """

rule pyclone_conn_post:
  input:
    loci="%s/result.loci.tsv" % OUTDIR,
    clust="%s/result.clusters.tsv" % OUTDIR
  output:
    snvs="%s/inf.snvs.csv" % OUTDIR,
    prev="%s/inf.clusters.csv" % OUTDIR
  shell:
    """
    bash {config[scripts]}/pyclone_parse_result.sh \
      {input.clust} {input.loci}
    """

rule pyclone_conn_metrics:
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
