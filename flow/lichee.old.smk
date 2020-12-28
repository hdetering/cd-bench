# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule lichee_prep:
  input:  "sim/snvs.filt.vcf"
  output: "sim/lichee/snv.vaf.tsv"
  params:
    dir = "lichee"
  log: "log/lichee_prep.log"
  shell:
    """
    module load miniconda3/4.7.10
    time (
    set -x
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/vcf2lichee.py \
      --vcf {input} \
      --out {output} \
      --normal "RN"
    ) >{log} 2>&1
    """

rule lichee:
  input:  "lichee/snv.vaf.tsv"
  output: "lichee/lichee.done"
  params:
    dir = "lichee"
  log: "log/lichee.log"
  shell:
    """
    module load jdk/8u181
    time (
    set -x
    {config[tools][lichee]} -build \
      -i {input} \
      --normal 3 \
      -maxVAFAbsent 0.01 \
      -minVAFPresent 0.05 \
      -maxVAFValid 1.0 \
      -dot \
      --verbose
    touch {output}
    ) >{log} 2>&1
    """

rule lichee_cn2_prep:
  input:  "sim/snv.cn2.vcf"
  output: "lichee/cn2/snv.vaf.tsv"
  params:
    dir = "lichee/cn2"
  log: "log/lichee_cn2_prep.log"
  shell:
    """
    time (
    set -x
    module load miniconda3/4.7.10
    if [[ ! -d "{params.dir}" ]]; then
      mkdir -p {params.dir}
    fi
    python {config[scripts]}/vcf2lichee.py \
      --vcf {input} \
      --out {output} \
      --normal "RN"
    ) >{log} 2>&1
    """

rule lichee_cn2:
  input:  "lichee/cn2/snv.vaf.tsv"
  output: "lichee/cn2/lichee.done"
  params:
    dir = "lichee/cn2"
  log: "log/lichee_cn2.log"
  shell:
    """
    module load jdk/8u181
    time (
    set -x
    {config[tools][lichee]} -build \
      -i {input} \
      --normal 0 \
      -maxVAFAbsent 0.01 \
      -minVAFPresent 0.05 \
      -maxVAFValid 0.6 \
      --verbose
    touch {output}
    ) >{log} 2>&1
    """

rule lichee_cn2_relaxed_prep:
  input:  "sim/snv.cn2.relaxed.vcf"
  output: "lichee/cn2_relaxed/snv.vaf.tsv"
  params:
    dir = "lichee/cn2_relaxed"
  log: "log/lichee_cn2_relaxed_prep.log"
  shell:
    """
    module load miniconda3/4.7.10
    time (
    set -x
    if [[ ! -d "{params.dir}" ]]; then
      mkdir -p {params.dir}
    fi
    python {config[scripts]}/vcf2lichee.py \
      --vcf {input} \
      --out {output} \
      --normal "RN"
    ) >{log} 2>&1
    """

rule lichee_cn2_relaxed:
  input:  "lichee/cn2_relaxed/snv.vaf.tsv"
  output: "lichee/cn2_relaxed/lichee.done"
  params:
    dir = "lichee/cn2_relaxed"
  log: "log/lichee_cn2_relaxed.log"
  shell:
    """
    module load jdk/8u181
    time (
    set -x
    {config[tools][lichee]} -build \
      -i {input} \
      --normal 0 \
      -maxVAFAbsent 0.01 \
      -minVAFPresent 0.05 \
      -maxVAFValid 0.6 \
      --verbose
    touch {output}
    ) >{log} 2>&1
    """

rule lichee_cn2_relaxed_post:
  input:  
    lichee = "lichee/cn2_relaxed/snv.vaf.tsv.trees.txt",
    #true_tree = "sim/true.tree.adj.csv",
    true_tree = "sim/clones.nex.adj.csv",
    true_snvs = "sim/true.snv_cluster.csv"
  output:
    trees = "lichee/cn2_relaxed/lichee.trees.csv",
    tree  = "lichee/cn2_relaxed/lichee.tree.0.csv",
    snvs  = "lichee/cn2_relaxed/lichee.snvs.csv",
    caset_u = "lichee/cn2_relaxed/lichee.tree.0.CASet_union.dist",
    caset_i = "lichee/cn2_relaxed/lichee.tree.0.CASet_isect.dist",
    disc_u = "lichee/cn2_relaxed/lichee.tree.0.DISC_union.dist",
    disc_i = "lichee/cn2_relaxed/lichee.tree.0.DISC_isect.dist"
  shell:
    """
    set -x
    python3 {config[scripts]}/lichee_parse_result.py {input.lichee}
    # extract first (and typically only) tree
    tail -n+2 {output.trees} | cut -f2,3 -d, > {output.tree}
    # calculate distance measures
    python3 {config[scripts]}/tree_dist.opt.py \
      --CASet --union --collapse \
      {input.true_tree} \
      {input.true_snvs} \
      {output.tree} \
      {output.snvs} \
    | tail -1 \
    | tee {output.caset_u}
    python3 {config[scripts]}/tree_dist.opt.py \
      --CASet --intersect --collapse \
      {input.true_tree} \
      {input.true_snvs} \
      {output.tree} \
      {output.snvs} \
    | tail -1 \
    | tee {output.caset_i}
    python3 {config[scripts]}/tree_dist.opt.py \
      --DISC --union --collapse \
      {input.true_tree} \
      {input.true_snvs} \
      {output.tree} \
      {output.snvs} \
    | tail -1 \
    | tee {output.disc_u}
    python3 {config[scripts]}/tree_dist.opt.py \
      --DISC --intersect --collapse \
      {input.true_tree} \
      {input.true_snvs} \
      {output.tree} \
      {output.snvs} \
    | tail -1 \
    | tee {output.disc_i}
    """
