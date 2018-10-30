mkdir 00_temp 0_fastq_raw 1_QC 2_fastq_trimmed 3_merged 4_aligned 5_flag 6_deduplicated 7_methCalls
mv *.fastq.gz 0_fastq_raw/
snakemake all
snakemake multiQC