# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule canopy_prep:
  input:
    csv="sim/snv.cn2.strict.csv",
  output:
    alt="canopy/rc_alt.csv",
    tot="canopy/rc_tot.csv"
  params:
    dir="canopy"
  log: "log/canopy_prep.log"
  group: "canopy"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/snvs2canopy.py \
      --csv {input.csv} \
      --normal "RN" \
      --outdir {params.dir}
    ) >{log} 2>&1
    """

rule canopy:
  input:
    alt="canopy/rc_alt.csv",
    tot="canopy/rc_tot.csv"
  output:
    snv="canopy/canopy.snvs.csv",
    tre="canopy/canopy.trees.csv",
    res="canopy/canopy.clusters.csv",
  log:
    "log/canopy.log"
  threads: 5
  group: "canopy"
  resources: mem_mb = 32000, time_mins = 240
  shell:
    """
    time (
   
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][canopy_env]}
    set -eu

    Rscript {config[scripts]}/run_canopy.R canopy {threads}

    ) >{log} 2>&1
    """

rule canopy_post:
  input:
    snv="canopy/canopy.snvs.csv",
    true="sim/true.snvs.csv"
  output:
    yml="canopy/metrics_clustering.yml"
  shell:
    """
    python {config[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {input.true} {input.snv} > {output.yml}
    """
