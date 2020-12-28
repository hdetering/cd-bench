# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

OUTPFX = "pairtree_cn2_seqz"

rule pairtree_prep:
  input:
    csv = "sequenza/snv.cn2.csv"
  output:
    snv = "%s/input.snvs.tsv" % OUTPFX,
    par = "%s/params.pre.json" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s_prep.log" % OUTPFX
  group: "pairtree"
  shell:
    """
    time (
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi
    python {config[scripts]}/snvs2pairtree.py \
      --csv {input.csv} \
      --outdir {params.dir}
    ) >{log} 2>&1
    """

rule pairtree_cluster:
  input:
    snv = "%s/input.snvs.tsv" % OUTPFX,
    par = "%s/params.pre.json" % OUTPFX
  output:
    par = "%s/params.json" % OUTPFX
  log: "log/%s_cluster.log" % OUTPFX
  threads: 4
  shell:
    """
    time (
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][pairtree_env]}
    set -eu

    set -x
    $CONDA_PREFIX/bin/python {config[tools][pairtree_dir]}/bin/clustervars \
      --parallel {threads} \
      {input.snv} \
      {input.par} \
      {output.par}
    ) >{log} 2>&1
    """

rule pairtree:
  input:
    snv = "%s/input.snvs.tsv" % OUTPFX,
    par = "%s/params.json" % OUTPFX
  output:
    done = "%s/inf.done" % OUTPFX
  params:
    out_npz = "%s/result.npz" % OUTPFX,
    outdir = OUTPFX
  log:
    "log/%s.log" % OUTPFX
  threads: 4
  group: "pairtree"
  shell:
    """
    time (
    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][pairtree_env]}
    set -eu

    set -x
    $CONDA_PREFIX/bin/python {config[tools][pairtree_dir]}/bin/pairtree \
      --parallel {threads} \
      --params {input.par} \
      {input.snv} \
      {params.out_npz}

    ) >{log} 2>&1

    touch {output.done}
    """

rule pairtree_post:
  input: "%s/inf.done" % OUTPFX
  output:
    snvs = "%s/inf.snvs.csv" % OUTPFX,
    prev = "%s/inf.clusters.csv" % OUTPFX,
    tree = "%s/inf.trees.csv" % OUTPFX
  params:
    in_snv = "%s/input.snvs.tsv" % OUTPFX,
    in_sum = "%s/result.summ.json.gz" % OUTPFX,
    #in_mut = "%s/result.muts.json.gz" % OUTPFX,
    in_ass = "%s/result.mutass.zip" % OUTPFX
  shell:
    """
    python3 {config[scripts]}/phylowgs_parse_result.py \
      {params.in_snv} {params.in_sum} {params.in_ass}
    """

rule pairtree_metrics:
  input:  
    inf_tree = "%s/inf.trees.csv" % OUTPFX,
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
