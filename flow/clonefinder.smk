# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
OUTPFX = "clonefinder_cn2_seqz"

rule clonefinder_prep:
  input:  "sequenza/snv.cn2.vcf"
  output: "%s/input.tsv" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s_prep.log" % OUTPFX
  shell:
    """
    module load miniconda3/4.7.10
    time (
    set -x
    if [[ ! -d "{params.dir}" ]]; then
      mkdir {params.dir}
    fi

    # copy required config files
    # (required when not cd'ing to CloneFinder dir for execution
    #cp {config[tools][clonefinder_dir]}/options.ini {params.dir}
    #cp {config[tools][clonefinder_dir]}/ancestral_seqs_MP_nucleotide.mao {params.dir}
    #cp {config[tools][clonefinder_dir]}/infer_MP_nucleotide.mao {params.dir}

    python {config[scripts]}/vcf2clonefinder.py \
      --vcf {input} \
      --out {output} \
      #--normal "RN"
    ) >{log} 2>&1
    """

rule clonefinder:
  input:  "%s/input.tsv" % OUTPFX
  output: 
    meg = "%s/inputsnv_CloneFinder.meg" % OUTPFX,
    nwk = "%s/inputsnv_CloneFinder.nwk" % OUTPFX,
    txt = "%s/inputsnv_CloneFinder.txt" % OUTPFX,
    log = "%s/inputsnv_summary.txt" % OUTPFX
  params:
    dir = OUTPFX
  log: "log/%s.log" % OUTPFX
  conda: "envs/clonefinder.yaml"
  shell:
    """
    export PATH=/mnt/netapp1/posadalab/APPS/mega/:$PATH
    time (
    set -x
    infile=$(readlink -f {input})
    cd {config[tools][clonefinder_dir]}
    which python
    python clonefinder.py \
      snv $infile
    ) >{log} 2>&1
    """

rule clonefinder_post:
  input:
    tsv = "%s/input.tsv" % OUTPFX,
    meg = "%s/inputsnv_CloneFinder.meg" % OUTPFX,
    nwk = "%s/inputsnv_CloneFinder.nwk" % OUTPFX,
    txt = "%s/inputsnv_CloneFinder.txt" % OUTPFX
  output:
    tree = "%s/inf.trees.csv" % OUTPFX,
    snvs = "%s/inf.snvs.csv" % OUTPFX,
    prev = "%s/inf.clusters.csv" % OUTPFX
  shell:
    """
    set -x

    python {config[scripts]}/clonefinder_parse_result.py\
      {input.tsv} {input.nwk} {input.meg} {input.txt}
    
    """

rule clonefinder_metrics:
  input:  
    inf_tree = "%s/inf.trees.csv" % OUTPFX,
    inf_snvs = "%s/inf.snvs.csv" % OUTPFX,
    inf_prev = "%s/inf.clusters.csv" % OUTPFX,
    true_tree = "sim/true.tree.csv",
    true_snvs = "sim/true.snvs.csv",
    true_prev = "sim/true.clusters.csv"
  output:
    clust = "%s/metrics_clustering.yml" % OUTPFX,
    prev  = "%s/metrics_prevalence.yml" % OUTPFX,
    phylo = "%s/metrics_phylogeny.yml" % OUTPFX
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

    # calculate tree distance measures
    # NOTE: add '--collapse' if running too slowly
    python3 {config[scripts]}/metrics_phylogeny.py \
      --CASet_isect --CASet_union --DISC_isect --DISC_union \
      --collapse --ignore-homoplasy \
      {input.true_tree} \
      {input.true_snvs} \
      {input.inf_tree} \
      {input.inf_snvs} \
    | tee {output.phylo}
    """
