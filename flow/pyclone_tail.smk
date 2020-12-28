# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

import os
import bisect
from glob import glob
import yaml
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.vcf.gz")
PWD = os.getcwd

rule pyclone_prep:
  input:
    vcf="sim/snv.cn2.tail.vcf",
  output:
    ["pyclone_tail/%s.tsv" % x for x in SAMPLES]
  params:
    dir="pyclone_tail"
  log: "log/pyclone_tail_prep.log"
  group: "pyclone"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/vcf2pyclone.py \
      --vcf {input.vcf} \
      --normal "RN" \
      --outdir {params.dir}
    ) >{log} 2>&1
    """

rule pyclone_yaml:
  input:  "pyclone_tail/{sample}.tsv"
  output: "pyclone_tail/{sample}.yaml"
  log:    "log/pyclone_tail_yaml.{sample}.log"
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

rule pyclone_conf:
  input:
    yml=["pyclone_tail/%s.yaml" % x for x in SAMPLES]
  output:
    cfg="pyclone_tail/config.yaml"
  group: "pyclone"
  run:
    import os, yaml
    cfg = {
      # Specifies working directory for analysis. All paths in the rest of the file are relative to this.
      'working_dir': os.path.join(os.getcwd(), 'pyclone_tail'),
      # Where the trace (output) from the PyClone MCMC analysis will be written.
      'trace_dir'  : 'trace',
      # Specifies which density will be used to model read counts. Most people will want pyclone_beta_binomial or pyclone_binomial
      'density'    : 'pyclone_beta_binomial',
      # Number of iterations of the MCMC chain.
      'num_iters'  : 10000,

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

rule pyclone:
  input:
    cfg="pyclone_tail/config.yaml"
  output:
    loci="pyclone_tail/result.loci.tsv",
    clust="pyclone_tail/result.clusters.tsv"
  log:
    "log/pyclone_tail.log"
  threads: 1
  group: "pyclone"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate pyclone
    set -eu
    PyClone run_analysis --config_file {input.cfg} 
    PyClone build_table --config_file {input.cfg} --out_file {output.loci} --table_type loci
    PyClone build_table --config_file {input.cfg} --out_file {output.clust} --table_type cluster

    ) >{log} 2>&1
    """

rule pyclone_post:
  input:
    loci="pyclone_tail/result.loci.tsv",
    clust="pyclone_tail/result.clusters.tsv",
    true="sim/true.snvs.csv"
  output:
    csv="pyclone_tail/pyclone.snvs.csv",
    yml="pyclone_tail/metrics_clustering.yml"
  shell:
    """
    bash {options[scripts]}/pyclone_parse_result.sh \
      {input.clust} {input.loci}

    python {options[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {input.true} {output.csv} > {output.yml}
    """
