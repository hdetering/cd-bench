# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule sequenza_vcf_hdr:
  input:  "sim/bam/RN.rc.vcf.gz"
  output: "vcf.header.txt"
  shell:
    """
    cat > {output}<< EOF
##FORMAT=<ID=CN1,Number=1,Type=Float,Description="True copy number of maternal allele.">
##FORMAT=<ID=CN2,Number=1,Type=Float,Description="True copy number of paternal allele.">
##FORMAT=<ID=CNI,Number=1,Type=Integer,Description="Inferred copy number.">
##FORMAT=<ID=CNIA,Number=1,Type=Integer,Description="Inferred copy number for allele A.">
##FORMAT=<ID=CNIB,Number=1,Type=Integer,Description="Inferred copy number for allele B.">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
EOF
    """

rule sequenza_gc:
  input:  "sim/ref.fa.gz"
  output: "sequenza/ref_gc50.wig.gz"
  log:    "log/sequenza_gc.log"
  threads: 1
  shell:
    """
    time (

    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][sequenza_env]}
    set -eu

    sequenza-utils gc_wiggle --fasta {input} -w 50 | gzip > {output}
    ) >{log} 2>&1
    """

rule sequenza_prep:
  input:
    vcf="sim/bam/{sample}.paired.vcf.gz",
    ref="sequenza/ref_gc50.wig.gz",
  output:
    srt="sequenza/{sample}.paired.sort.vcf.gz",
    snp="sequenza/{sample}.snps.vcf",
    seqz="sequenza/{sample}.seqz.gz"
  log:  "log/{sample}.sequenza_prep.log"
  threads: 1
  shell:
    """
    time(
    module load bcftools
    bcftools sort -Oz {input.vcf} > {output.srt}

    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][sequenza_env]}
    set -eu

    zcat {input.vcf} | awk -f {config[scripts]}/vcf_fix_AD.awk \
    | awk '/^#/||($3~/^g/&&$3!~/,/&&$10!~/^\./&&$11!~/^\./)' \
    > {output.snp}

    sequenza-utils snp2seqz --vcf {output.snp} -gc {input.ref} \
    | gzip -c > {output.seqz}
    
    ) >{log} 2>&1
    """

rule sequenza:
  input:  "sequenza/{sample}.seqz.gz"
  output: "sequenza/{sample}_segments.txt"
  log:    "log/{sample}.sequenza.log"
  shell:
    """
    time(

    # to fix https://github.com/conda/conda/issues/8186
    set +eu
    source activate {config[tools][sequenza_env]}
    set -eu

    Rscript {config[scripts]}/run_sequenza.R {input}

    ) >{log} 2>&1
    """

rule sequenza_post:
  input:
    hdr="vcf.header.txt",
    vcf="sim/bam/{sample}.filt.cn.vcf.gz",
    seg="sequenza/{sample}_segments.txt"
  output:
    tsv="sequenza/{sample}.anno.tsv.gz",
    vcf="sequenza/{sample}.cn.inf.vcf.gz"
  log:  "log/{sample}.sequenza_post.log"
  shell:
    """
    module load bedtools bcftools

    time(
    
    bedtools intersect -a {input.vcf} -b <(tail -n+2 {input.seg}) -wb \
    | cut -f 1,2,21,22 \
    | bgzip > {output.tsv}
    tabix -s1 -b2 -e2 -c"#" {output.tsv}

    bcftools annotate -a {output.tsv} -c CHROM,POS,FMT/CNIA,FMT/CNIB -h {input.hdr} -Oz {input.vcf} \
    > {output.vcf}
    bcftools index {output.vcf}

    ) >{log} 2>&1
    """

# Merge VCF files from multiple samples into single VCF
rule sequenza_merge:
  input:
    vcf=["sequenza/{}.cn.inf.vcf.gz".format(s) for s in SAMPLES_TUM]
  output: 
    vcf = "sequenza/snv.cn.vcf",
    csv = "sequenza/snv.cn.csv"
  params:
    var=" ".join(["sequenza/{}.cn.inf.vcf.gz".format(s) for s in SAMPLES_TUM])
  log:
    "log/sequenza_merge.log"
  shell:
    """
    module load bcftools
    time (
    set -x
    bcftools merge --force-samples {params.var} \
    | awk '/^#/{{print;next}}$3!~/^g/{{print}}' \
    | awk -f {config[scripts]}/vcf_fix_AD.awk \
    > {output.vcf}

    python3 {config[scripts]}/vcf2csv.seqz.py {output.vcf} > {output.csv}
    ) 1>{log} 2>&1
    """

# Filter multi-sample for variants with all samples having inferred copy number =2
rule sequenza_filt_cn2:
  input:  "sequenza/snv.cn.vcf"
  output: 
    vcf="sequenza/snv.cn2.vcf",
    csv="sequenza/snv.cn2.csv"
  shell:
    """
    module load bcftools
    bcftools view -e 'FMT/CNI!=2 & FMT/CNI!="."' {input} > {output.vcf}

    python3 {config[scripts]}/vcf2csv.py {output.vcf} > {output.csv}
    """
