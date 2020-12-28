# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule lichee_tail_prep:
  input:  "sim/snv.cn2.tail.vcf"
  output: "lichee_tail/snv.vaf.tsv"
  params:
    dir = "lichee_tail"
  log: "log/lichee_tail_prep.log"
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

rule lichee_tail:
  input:  "lichee_tail/snv.vaf.tsv"
  output: "lichee_tail/lichee_tail.done"
  params:
    dir = "lichee_tail"
  log: "log/lichee_tail.log"
  shell:
    """
    module load jdk/8u181
    time (
    set -x
    {config[tools][lichee]} -build \
      -i {input} \
      --normal 0 \
      -maxVAFAbsent 0.05 \
      -minVAFPresent 0.10 \
      -maxVAFValid 0.7 \
      --verbose
    touch {output}
    ) >{log} 2>&1
    """

rule lichee_tail_post:
  input:  
    lichee = "lichee_tail/snv.vaf.tsv.trees.txt"
  output:
    trees = "lichee_tail/lichee_tail.trees.csv",
    tree  = "lichee_tail/lichee_tail.tree.0.csv",
    snvs  = "lichee_tail/lichee_tail.snvs.csv"
  shell:
    """
    set -x
    python3 {config[scripts]}/lichee_parse_result.py {input.lichee}
    # extract first (and typically only) tree
    tail -n+2 {output.trees} | cut -f2,3 -d, > {output.tree}
    """

rule lichee_tail_metr_clust:
  input:  
    tree = "lichee_tail/lichee_tail.tree.0.csv",
    snvs = "lichee_tail/lichee_tail.snvs.csv",
    true_tree = "sim/clones.nex.adj.csv",
    true_snvs = "sim/true.snvs.csv"
  output:
    "lichee_tail/metrics_phylo.yml"
  shell:
    """
    set -x
    # calculate distance measures
    # NOTE: add '--collapse' if running too slowly
    python3 {config[scripts]}/tree_dist.opt.py \
      --CASet_isect --CASet_union --DISC_isect --DISC_union --collapse\
      {input.true_tree} \
      {input.true_snvs} \
      {input.tree} \
      {input.snvs} \
    | tee {output}
    """
