# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule spruce_prep:
  input:  "sim/snv.cn2.strict.vcf"
  output: "spruce/input.tsv"
  params:
    dir = "spruce"
  log: "log/spruce_prep.log"
  shell:
    """
    module purge
    module load miniconda3/4.7.10
    # fix to avoid error
    # https://github.com/conda/conda/issues/8186#issuecomment-532874667
    set +eu
    source activate {config[tools][spruce_env]}
    set -eu

    time (
    set -x
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi

    python {config[scripts]}/vcf2spruce.py \
      --vcf {input} \
      --out {output} \
      --normal "RN"
    ) >{log} 2>&1
    """

rule spruce_cliques:
  input:  "spruce/input.tsv"
  output: "spruce/input.cliques"
  log: "log/spruce_cliques.log"
  shell:
    """
    module purge
    module load lemon/1.3.1
    module load gcc/6.4.0 openmpi/2.1.1 boost/1.68.0-python-2.7.15
    
    time (
    set -x
    {config[tools][spruce_dir]}/cliques {input} > {output}
    ) >{log} 2>&1
    """

rule spruce_enumerate:
  input:
    tsv = "spruce/input.tsv",
    clq = "spruce/input.cliques"
  output: 
    res = "spruce/input.res.gz",
    mrg = "spruce/input.merged.res"
  log: "log/spruce_enumerate.log"
  params:
    time_secs = 28800
  threads: 4
  shell:
    """
    module purge
    module load lemon/1.3.1
    module load gcc/6.4.0 openmpi/2.1.1 boost/1.68.0-python-2.7.15
    
    time (
    set -x
    {config[tools][spruce_dir]}/enumerate \
      {input.tsv} -clique {input.clq} \
      -t {threads} -ll {params.time_secs} \
    | gzip > {output.res}

    gzcat {output.res} \
    | {config[tools][spruce_dir]}/rank - \
    > {output.mrg}
    ) >{log} 2>&1
    """
rule spruce_post:
  input:
    meg = "spruce/inputsnv_CloneFinder.meg",
    nwk = "spruce/inputsnv_CloneFinder.nwk",
    txt = "spruce/inputsnv_CloneFinder.txt"
  output:
    adj = "spruce/spruce.tree.adj.csv",
    snv = "spruce/spruce.snv_cluster.csv",
    pre = "spruce/spruce.cluster_prev.csv"
  shell:
    """
    set -x

    python {config[scripts]}/tree_to_ajacency.py\
      --newick {input.nwk} \
    > {output.adj}

    
    """
