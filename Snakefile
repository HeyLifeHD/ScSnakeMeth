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
        

rule multiQC:
    message:
        "Perform MultiQC"
    shell:
        """
        multiqc --title "Raw FastQC reports" --filename Raw_FastQC_report --outdir 1_MultiQC/ ./
        """ 


rule trim:
    message:
        "Preprocessing of fastq files"
    input:
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
        bam="4_aligned/{sample}_merged_bismark_bt2.bam",
        log="4_aligned/{sample}_bismark_bt2.bam.log"
    params: 
        ref= BISMARK_REFERENCE
    shell:
        """bismark --non_directional --parallel 2 --output_dir 4_aligned/ \
        --rg_id {sample} --rg_sample {sample} -p 4 {params.ref} {input} > {output.log}
        touch {output.bam}
        """


rule sort:
    message:
        "Sorting bam files"
    input:
        "4_aligned/{sample}_merged_bismark_bt2.bam"
    output:
        "5_flag/{sample}_sorted.bam"
    shell:
        "samtools sort -@ 4 -o {output} {input}"


rule dedup2:
    message:
        "Remove duplicates"
    input:
        "5_flag/{sample}_sorted.bam"
    output:
        idx="5_flag/{sample}_sorted.bam.idxstats",
        bam="6_deduplicated/{sample}_deduplicated.bam",
        metric="6_deduplicated/{sample}_deduplicated.bam.picard.metric"
    shell:
        """
        samtools index {input}
        samtools idxstats  {input} > {output.idx} 
        picard MarkDuplicates I= {input} O= {output.bam} M={output.metric} CREATE_INDEX= true
    
        """


rule methExtract:
    message:
        "Extract methylation calls"
    input:
        "6_deduplicated/{sample}_deduplicated.bam"
    output:
        "7_methCalls/{sample}_bismark_bt2.deduplicated.bedGraph.gz"
    shell:
        """
        bismark_methylation_extractor -s --comprehensive --merge_non_CpG -o 7_methCalls/ --parallel 6 --bedGraph --gzip {input}
        touch {output}
        """         
