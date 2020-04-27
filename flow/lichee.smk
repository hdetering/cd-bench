# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

SAMPLES, = glob_wildcards("sim/bam/{sample,R\d+}.rc.filt.vcf")

rule lichee_prep:
  input:
  output:
  params:
    dir = "lichee"
  shell:
    """
    time (
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    ) >{log} 2>&1
    """

rule lichee:
