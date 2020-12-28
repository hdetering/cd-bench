# vim: syntax=python tabstop=2 expandtab shiftwidth=2
# coding: utf-8

# VCF header definitions for FORMAT and INFO tags
# rule vcf_header:
#   input:  "sim/bam/RN.rc.vcf.gz"
#   output: "vcf.header.txt"
#   shell:
#     """
#     cat > {output}<< EOF
# ##FORMAT=<ID=CN1,Number=1,Type=Float,Description="True copy number of maternal allele.">
# ##FORMAT=<ID=CN2,Number=1,Type=Float,Description="True copy number of paternal allele.">
# ##FORMAT=<ID=CNI,Number=1,Type=Integer,Description="Inferred copy number.">
# ##FORMAT=<ID=CNIA,Number=1,Type=Integer,Description="Inferred copy number for allele A.">
# ##FORMAT=<ID=CNIB,Number=1,Type=Integer,Description="Inferred copy number for allele B.">
# ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
# EOF
#     """

# annotate VCF with true copy number
rule vcf_anno_cn_true:
  input:
    vcf = "sim/bam/{sample}.rc.filt.vcf.gz",
    bed = "sim/bed/{sample}.cn.bed.gz",
    hdr = "vcf.header.txt"
  output:
    vcf = "sim/bam/{sample}.filt.cn.vcf.gz",
    idx = "sim/bam/{sample}.filt.cn.vcf.gz.csi"
  log: "log/{sample}.vcf_anno_cn_true.log"
  shell:
    """
    module load bcftools
    time (
    set -x
    bcftools annotate \
      -a {input.bed} \
      -c CHROM,FROM,TO,FMT/CN1,FMT/CN2 \
      -h {input.hdr} \
      -Oz {input.vcf} \
    > {output.vcf}
    bcftools index {output.vcf}
    ) >{log} 2>&1
    """

rule vcf_merge:
  input:
    tum = "sim/bam/{sample}.filt.cn.vcf.gz",
    nrm = "sim/bam/RN.filt.cn.vcf.gz"
  output:
    vcf = "sim/bam/{sample}.paired.vcf.gz",
    idx = "sim/bam/{sample}.paired.vcf.gz.csi"
  log: "log/{sample}.vcf_merge.log"
  shell:
    """
    module load bcftools
    time (
    set -x
    bcftools merge {input.nrm} {input.tum} -Oz > {output.vcf} &&\
    bcftools index {output.vcf}
    ) >{log} 2>&1
    """

# calculate BAF and LRR for VCF
rule vcf_calc_baf_lrr:
  input:  "sim/bam/{sample}.paired.vcf.gz"
  output: "sim/bam/{sample}.baf.vcf.gz"
  log: "log/{sample}.vcf_calc_baf_lrr.log"
  shell:
    """
    module load bcftools
    time(
    set -x
    $CONDA_PREFIX/bin/python {config[scripts]}/vcf_calc_baf_lrr.py {input} RN \
    | bgzip -c > {output} &&\
    bcftools index {output}
    ) >{log} 2>&1
    """

# infer CN state
rule bcftools_cnv:
  input:  "sim/bam/{sample}.baf.vcf.gz"
  output:
    tab = "bcftools_cnv/{sample}/cn.paired.tab.gz",
    idx = "bcftools_cnv/{sample}/cn.paired.tab.gz.tbi"
  params:
    outdir = "bcftools_cnv/{sample}"
  log: "log/{sample}.bcftools_cnv.log"
  shell:
    """
    module load bcftools
    time(
    set -x

    bcftools cnv -c RN -s {wildcards.sample} -p 0 -o {params.outdir} {input}
    
    # post-process output (create annotation file)
    awk -v OFS="\t" 'NR==FNR{{m[$1" "$2]=$3;next}}{{print $1,$2,m[$1" "$2],$3}}' \
      {params.outdir}/cn.RN.tab {params.outdir}/cn.{wildcards.sample}.tab \
    | bgzip -c > {output.tab}
    tabix -s1 -b2 -e2 -c"#" {output.tab}

    ) >{log} 2>&1
    """

# annotate VCF with inferred copy number
rule vcf_anno_cn_inf:
  input:  
    vcf = "sim/bam/{sample}.baf.vcf.gz",
    tab = "bcftools_cnv/{sample}/cn.paired.tab.gz"
  output:
    vcf = "bcftools_cnv/{sample}.cn.inf.vcf.gz",
    idx = "bcftools_cnv/{sample}.cn.inf.vcf.gz.csi"
  log: "log/{sample}.vcf_anno_cn_inf.log"
  shell:
    """
    module load bcftools
    time(
    set -x
    bcftools annotate -a {input.tab} -c CHROM,POS,FMT/CNI -Oz {input.vcf} > {output.vcf} &&\
    bcftools index {output.vcf}
    ) >{log} 2>&1
    """

# Merge VCF files from multiple samples into single VCF
rule bcftools_merge:
  input:
    vcf=["bcftools_cnv/{}.cn.inf.vcf.gz".format(s) for s in SAMPLES_TUM]
  output: 
    "bcftools_cnv/snv.cn.vcf"
  params:
    var=" ".join(["bcftools_cnv/{}.cn.inf.vcf.gz".format(s) for s in SAMPLES_TUM])
  log:
    "log/combine_variants.log"
  shell:
    """
    module load bcftools
    time (
    set -x
    bcftools merge --force-samples {params.var} \
    | awk '/^##/{{print;next}}/^#CHROM/||$3!~/^g/{{printf $1;for(i=2;i<=NF;i++){{if(i<=11||i%2==1)printf "\\t"$i}};printf "\\n"}}' \
    | awk -f {config[scripts]}/vcf_fix_AD.awk \
    > {output}
    ) 1>{log} 2>&1
    """

# Filter multi-sample for variants with all samples having inferred copy number =2
rule bcftools_filt_cn2:
  input:  "bcftools_cnv/snv.cn.vcf"
  output: 
    vcf = "bcftools_cnv/snv.cn2.vcf",
    csv = "bcftools_cnv/snv.cn2.csv"
  shell:
    """
    module load bcftools
    bcftools view -e 'FMT/CNI!=2 & FMT/CNI!="."' {input} > {output.vcf}

    python3 {config[scripts]}/vcf2csv.py {output.vcf} > {output.csv}
    """
