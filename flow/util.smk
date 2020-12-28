# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

# index FASTA file
rule ref_index:
  input:
    "{ref}.fa"
  output:
    fai="{ref}.fa.fai",
    dict="{ref}.dict"
  log:    "log/{ref}.ref_index.log"
  threads: 1
  shell:
    """
    time (

    module load gcc/6.4.0 
    module load samtools/1.9
    module load jdk/8u181

    samtools faidx {input}
    {config[tools][picard]} CreateSequenceDictionary \
      REFERENCE={input} \
      OUTPUT={output.dict}

    ) >{log} 2>&1
    """

# Merge VCF files from multiple samples into single VCF
rule combine_cn_vcfs:
  input:
    nrm="bcftools_cnv/RN.cn.inf.vcf.gz"
  output: 
    "sim/snvs.cn.vcf.gz"
  log:
    "log/combine_cn_vcfs.log"
  shell:
    """
    module load bcftools
    time (
    fn=$(basename {input.nrm})
    python {config[scripts]}/combine_vcfs.py $(
      for f in ${{{input.nrm}//RN/R*}}; do
        fn=$(basename $f);
        printf " --vcf %s:%s" ${fn%.cn.inf.vcf.gz} $f;
      done)
      --out {output}
    ) 1>{log} 2>&1
    """
#rule combine_cn_vcfs:
#  input:
#    nrm="bcftools_cnv/RN.cn.inf.vcf.gz"
#  output: 
#    "sim/snvs.cn.vcf.gz"
#  log:
#    "log/combine_cn_vcfs.log"
#  shell:
#    """
#    time (
#    fn=$(basename {input.nrm})
#    python {config[scripts]}/combine_vcfs.py $(
#      for f in ${{{input.nrm}//RN/R*}}; do
#        fn=$(basename $f);
#        printf " --vcf %s:%s" ${fn%.cn.inf.vcf.gz} $f;
#      done)
#      --out {output}
#    ) 1>{log} 2>&1
#    """

# Create sequence dictionary for reference
rule gatk_create_seq_dict:
  input:  "sim/ref.fa.gz"
  output: "sim/ref.dict"
  shell:
    """
    module load gatk/4.1.6.0
    gatk CreateSequenceDictionary -R {input}
    """

# Change contig records in VCF header
rule vcf_header_replace_contigs:
  input:  "bcftools_cnv/{sample}.cn.inf.vcf.gz"
  output: "bcftools_cnv/{sample}.cn.inf.mod.vcf"
  shell:
    """
    fn={input}
    zcat $fn > ${{fn%.gz}}
    awk -f {config[scripts]}/vcf_header_replace_contigs.awk \
      {config[resources]}/vcf.header.contig.txt \
      ${{fn%.gz}} \
    > {output}
    """

# Sort VCF files
rule bcftools_sort:
  input:  "bcftools_cnv/{sample}.cn.inf.mod.vcf"
  output: "bcftools_cnv/{sample}.cn.inf.sort.vcf"
  shell:
    """
    module load bcftools
    bcftools sort {input} > {output}
    """

# Merge VCF files from multiple samples into single VCF
rule combine_variants_tail:
  input:
    vcf=["bcftools_cnv/{}.cn.inf.vcf.gz".format(s) for s in SAMPLES_TUM]
  output: 
    "sim/snv.cn.tail.vcf"
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
# Merge VCF files from multiple samples into single VCF
#rule combine_variants:
#  input:
#    ref="sim/ref.fa",
#    vcf=["bcftools_cnv/{}.cn.inf.sort.vcf".format(s) for s in SAMPLES_TUM]
#  output: 
#    "sim/snv.cn.vcf"
#  params:
#    var=" ".join(["-V bcftools_cnv/{}.cn.inf.sort.vcf".format(s) for s in SAMPLES_TUM])
#  log:
#    "log/combine_variants.log"
#  shell:
#    """
#    module load jdk/8u181
#    time (
#      {config[tools][gatk3]} -T CombineVariants \
#        -R {input.ref} \
#        {params.var} \
#        -o {output} \
#        -genotypeMergeOptions UNSORTED
#    ) 1>{log} 2>&1
#    """

# Filter multi-sample for variants with all samples having inferred copy number =2
rule vcf_filt_cn2_strict:
  input:  "sim/snv.cn.vcf"
  output: "sim/snv.cn2.strict.vcf"
  shell:
    """
    module load bcftools
    bcftools view -e 'FMT/CNI!=2 & FMT/CNI!="."' {input} > {output}
    """ 

# Filter multi-sample for variants, exclude samples having inferred copy number !=2
rule vcf_filt_cn2_relaxed:
  input:  "sim/snv.cn.vcf"
  output: "sim/snv.cn2.relaxed.vcf"
  shell:
    """
    awk -f {config[scripts]}/vcf_filt_cn2_relaxed.awk {input} > {output}
    """ 

# create pileup from bam
rule bam2pileup:
    input:
        bam="{sample}.bam",
        ref="ref.fa"
    output: "{sample}.pu.gz"
    log:    "log/{sample}.bam2pileup.log"
    threads: 1
    shell:
        """
        module purge
        module load bcftools/1.9

        time (
        bcftools mpileup -f {input.ref} -Q 20 {input.bam} \
        | gzip > {output}
        ) >{log} 2>&1
        """

# Add read group to BAM files (required by GATK tools)
rule add_rg:
    input:  "data/{rep}/{sam}.bam"
    output: "data/{rep}/{sam}.RG.bam"
    log:    "data/{rep}/log/{sam}.picard_AddRG.log"
    group:  "mutect1"
    shell:
        """
        module load jdk/9.0.4
        java -jar {config[tools][picard]} AddOrReplaceReadGroups \
            I={input} \
            OUTPUT={output} \
            RGID={wildcards.sam} \
            RGSM={wildcards.sam} \
            RGPL="ILLUMINA" \
            RGLB="ART_500_pe" \
            RGPU="simulation" \
            CREATE_INDEX=true
        """

# Remove PCR duplicates
rule rm_dup:
    input:  "data/{rep}/{sam}.bam"
    output: "data/{rep}/{sam}.rmdup.bam"
    log:    "data/{rep}/log/{sam}.picard_MarkDup.log"
    params:
        tmp_dir="/tmp/{rep}",
        dupfile="data/{rep}/{sam}.dup.txt"
    shell:
        """
        module load jdk/9.0.4
        #module load picard/2.2.1

        time(
        java -jar {config[tools][picard]} MarkDuplicates \
            INPUT={input} \
            OUTPUT={output} \
            CREATE_INDEX=true \
            REMOVE_DUPLICATES=true \
            TMP_DIR={params.tmp_dir} \
            M={params.dupfile} \
            VALIDATION_STRINGENCY=LENIENT \
        > {log} 2>&1
        )
        """

rule bedtools_coverage:
  input:
      bam="{sample}.bam"
  output:
    "bedtools/{sample}.bam.cvg.bed"
  shell:
    """
    """    
