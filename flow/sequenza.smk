# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule sequenza_gc:
  input:  "ref.fa"
  output: "ref.gc50Base.txt.gz"
  log:    "log/sequenza.gc50.log"
  threads: 1
  shell:
    """
    time (
    config[tools][sequenza-utils] GC-windows -w 50 {input} | gzip > {output}
    ) >{log} 2>&1
    """
  

rule sequenza_prep:
  input:
    tumour="{sample}.pu.gz",
    normal="RN.pu.gz",
    ref="ref.fa",
    ref_gc50="ref.gc50Base.txt.gz"
  output:
    tumour="{sample}.seqz"
  log:  "log/sequenza.prep.log"
  threads: 1
  shell:
    """
    time() >{log} 2>&1
    """
