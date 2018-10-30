import sys

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")


import glob
import re

BISMARK_REFERENCE = "/bigdisk/ref_db/mm10/bismark/"
SAMPLES, = glob_wildcards("0_fastq_raw/{sam}_R1.fastq.gz")


for sample in SAMPLES:
  message("sample" + str(sample) + " will be processed")
 

rule all:
    input:
        expand("7_methCalls/{sample}/{sample}_bismark_bt2.deduplicated.bedGraph.gz", sample=SAMPLES)
#        expand("4_aligned/{sample}/{sample}_bismark_bt2.bam", sample=SAMPLES)
#        expand("4_aligned/{sample}/", sample=SAMPLES)
#        expand("5_flag/{sample}/{sample}_sorted.bam", sample=SAMPLES)
#        expand("4_aligned/{sample}/{sample}.log", sample=SAMPLES)
#        expand("3_merged/{sample}_merged.fastq.gz", sample=SAMPLES)

#rule dirCreate:
#    message:
#        "Create directories"
#    shell:
#        "mkdir 00_temp 1_QC 2_fastq_trimmed 3_merged 4_aligned 5_flag 6_deduplicated 7_methCalls"

#rule dirCreate:
#    message:
#        "Create directories"
#    shell:
#        "mkdir 00_temp 0_fastq_raw 1_QC 2_fastq_trimmed 3_merged 4_aligned 5_flag 6_deduplicated 7_methCalls"


#rule move:
#    message:
#        "Fastq files are moved into 0_fastq_raw"
#    input:
#        fwd="{sample}_R1.fastq.gz",
#        rev="{sample}_R2.fastq.gz"
#        "*.fastq.gz"
#    output:
#        fwd="0_fastq_raw/{sample}_R1.fastq.gz",
#        rev="0_fastq_raw/{sample}_R2.fastq.gz"
#        "0_fastq_raw/{*.fastq.gz"
#    shell:
#       """
#       #mv {input} {output}
#       mv {input.fwd} {output.fwd}
#       mv {input.rev} {output.rev}
#       """


rule multiQC:
    message:
        "Perform MultiQC"
    shell:
        "multiqc --title Raw FastQC reports --filename Raw_FastQC_report --outdir 1_MultiQC/" 


rule trim:
    message:
        "Preprocessing of fastq files"
    input:
#        fwd="{sample}_R1.fastq.gz",
#        rev="{sample}_R2.fastq.gz"
        fwd="0_fastq_raw/{sample}_R1.fastq.gz",
        rev="0_fastq_raw/{sample}_R2.fastq.gz"
    output:
        fwd="2_fastq_trimmed/{sample}_R1_trimmed.fastq.gz",
        rev="2_fastq_trimmed/{sample}_R2_trimmed.fastq.gz"
    shell:
        """
        fastp -i {input.fwd} -o {output.fwd} -I {input.rev} -O {output.rev} --trim_front1 10 --trim_front2 10 \
        --low_complexity_filter --thread 4 --verbose --html 1_QC/{sample}_fastp.html -R {sample} \
        --json 1_QC/{sample}_fastp.json -R {sample}
        """
        #mv {input.fwd} 0_fastq_raw/{sample}_R1.fastq.gz"
        #mv {input.rev} 0_fastq_raw/{sample}_R2.fastq.gz


rule merge:
    message:
        "Merging of fastq files"
    input:
        fwd="2_fastq_trimmed/{sample}_R1_trimmed.fastq.gz",
        rev="2_fastq_trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        "3_merged/{sample}_merged.fastq.gz"
    shell:
        "cat {input.fwd} {input.rev} > {output}"


rule align:
    message:
        "Aligning of fastaq files"
    input: 
        "3_merged/{sample}_merged.fastq.gz"
    output:
        "4_aligned/{sample}/{sample}_bismark_bt2.bam"
    params: 
        ref= BISMARK_REFERENCE
    shell:
        """bismark --non_directional --parallel 2 --output_dir 4_aligned/{sample}/ \
        --rg_id {sample} --rg_sample {sample} -p 4 {params.ref} {input} 
        touch {output}
        """


rule sort:
    message:
        "Sorting bam files"
    input:
        "4_aligned/{sample}/{sample}_bismark_bt2.bam"
    output:
        "5_flag/{sample}/{sample}_sorted.bam"
    shell:
        "samtools sort -@ 4 -o {output} {input}"


rule flagDup:
    message:
        "Flag and index duplicated reads"
    input:
        "5_flag/{sample}/{sample}_sorted.bam"
    output:
        "5_flag/{sample}/{sample}_sorted.bam.flagstat"
    shell:
        "samtools flagstat -@ 4 {input} > {output}"

#rule flagDup:
#    message:
#        "Flag and index duplicated reads"
#    input:
#        "5_flag/{sample}/{sample}_sorted.bam"
#    output:
#        "00_temp/{sample}_sorted.bam.flagstat"
#    shell:
#        "samtools flagstat -@ 4 {input} > {output}"

rule idx:
    message:
        "Flag and index duplicated reads" 
    input:
        "5_flag/{sample}/{sample}_sorted.bam"
    output:
        "5_flag/{sample}/{sample}_sorted.bam.idxstats"
    shell:
        "samtools idxstats {input} > {output}"

#rule idx:
#    message:
#        "Flag and index duplicated reads" 
#    input:
#        "5_flag/{sample}/{sample}_sorted.bam"
#    output:
#        "00_temp/{sample}_sorted.bam.idxstats"
#    shell:
#        "samtools idxstats {input} > {output}"

#rule cpidx:
#    message:
#        "Copy index files" 
#    input:
#        "00_temp/{sample}_sorted.bam.idxstats"
#    output:
#        "5_flag/{sample}/{sample}_sorted.bam.idxstats"
#    shell:
#        "cp {input} > {output}"      

rule dedup:
    message:
        "Remove duplicates"
    input:
        "5_flag/{sample}/{sample}_sorted.bam"
    output:
        "6_deduplicated/{sample}_deduplicated.bam"
    shell:
        "deduplicate_bismark --bam --output_dir {output} -s {input}"


rule methExtract:
    message:
        "Extract methylation calls"
    input:
        "6_deduplicated/{sample}_deduplicated.bam"
    output:
        "7_methCalls/{sample}/{sample}_bismark_bt2.deduplicated.bedGraph.gz"
    shell:
        """bismark_methylation_extractor -s --comprehensive --merge_non_CpG -o 7_methCalls/{sample}/ --parallel 6 --bedGraph --gzip {input}
        touch {output}"""
