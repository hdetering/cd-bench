# vim: syntax=python tabstop=2 expandtab shiftwidth=2
# coding: utf-8

# VCF header definitions for FORMAT and INFO tags
rule vcf_header:
  input:  "sim/bam/RN.rc.vcf.gz"
  output: "vcf.header.txt"
  shell:
    """
    cat > {output}<< EOF
##FORMAT=<ID=CN1,Number=1,Type=Float,Description="True copy number of maternal allele.">
##FORMAT=<ID=CN2,Number=1,Type=Float,Description="True copy number of paternal allele.">
##FORMAT=<ID=CNI,Number=1,Type=Integer,Description="Inferred copy number.">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
EOF
    """

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
    python {config[scripts]}/vcf_calc_baf_lrr.py {input} RN \
    | bgzip -c > {output} &&\
    bcftools index {output}
    ) >{log} 2>&1
    """

# infer CN state
rule bcftools_cnv_tail:
  input:  "sim/bam/{sample}_tail.baf.vcf.gz"
  output:
    tab = "bcftools_cnv/{sample}_tail/cn.paired.tab.gz",
    idx = "bcftools_cnv/{sample}_tail/cn.paired.tab.gz.tbi"
  params:
    outdir = "bcftools_cnv/{sample}_tail"
  log: "log/{sample}_tail.bcftools_cnv.log"
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
