# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
OUTPFX = "lichee_cn2_seqz"

rule lichee_prep:
  input:  "sequenza/snv.cn2.vcf"
  output: "%s/snv.vaf.tsv" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s_prep.log" % OUTPFX
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
    #  --normal "RN"
    ) >{log} 2>&1
    """

rule lichee:
  input:  "%s/snv.vaf.tsv" % OUTPFX
  output: "%s/inf.done" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s.log" % OUTPFX
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

rule lichee_post:
  input:  
    lichee = "%s/snv.vaf.tsv.trees.txt" % OUTPFX
  output:
    trees = "%s/inf.trees.csv" % OUTPFX,
    tree  = "%s/inf.tree.0.csv" % OUTPFX,
    snvs  = "%s/inf.snvs.csv" % OUTPFX,
    prev  = "%s/inf.clusters.csv" % OUTPFX
  shell:
    """
    set -x
    python3 {config[scripts]}/lichee_parse_result.py {input.lichee}
    # extract first (and typically only) tree
    tail -n+2 {output.trees} | cut -f2,3 -d, > {output.tree}
    """

rule lichee_metrics:
  input:  
    inf_tree = "%s/inf.tree.0.csv" % OUTPFX,
    inf_snvs = "%s/inf.snvs.csv" % OUTPFX,
    inf_prev = "%s/inf.clusters.csv" % OUTPFX,
    true_tree = "sim/true.tree.csv",
    true_snvs = "sim/true.snvs.csv",
    true_prev = "sim/true.clusters.csv"
  output:
    clust = "%s/metrics_clustering.yml" % OUTPFX,
    tree  = "%s/metrics_phylogeny.yml" % OUTPFX,
    prev  = "%s/metrics_prevalence.yml" % OUTPFX
  shell:
    """
    set -x
    python3 {config[scripts]}/metrics_clustering.py \
      --ARI --V-measure \
      {input.true_snvs} {input.inf_snvs} \
    | tee {output.clust}

    # calculate distance measures
    # NOTE: add '--collapse' if running too slowly
    python3 {config[scripts]}/metrics_phylogeny.py \
      --CASet_isect --CASet_union --DISC_isect --DISC_union --collapse\
      {input.true_tree} \
      {input.true_snvs} \
      {input.inf_tree} \
      {input.inf_snvs} \
    | tee {output.tree}

    python3 {config[scripts]}/metrics_prevalence.py \
      {input.true_prev} \
      {input.true_snvs} \
      {input.inf_prev} \
      {input.inf_snvs} \
    | tee {output.prev}
    """
