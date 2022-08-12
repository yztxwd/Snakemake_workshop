rule all:
    input:
        "output/mapped/final.merge.bam",
        "output/multiqc/multiqc_report.html"

#fastqc1
rule fastqc:
    input:
        "data/samples/SRR922442_1.fastq"
    output:
        html="output/fastqc/SRR922442_1_fastqc.html",
        zip="output/fastqc/SRR922442_1_fastqc.zip"
    shell:
        "fastqc {input} -o output/fastqc"

#fastqc2
#rule fastqc_wildcards:
#    input:
#        "data/samples/{prefix}.fastq"
#    output:
#        html="output/fastqc/{prefix}_fastqc.html",
#        zip="output/fastqc/{prefix}_fastqc.zip"
#    shell:
#        "fastqc {input} -o output/fastqc"

#fastqc3
rule fastqc_log:
    input:
        "data/samples/{prefix}.fastq"
    output:
        html="output/fastqc/{prefix}_fastqc.html",
        zip="output/fastqc/{prefix}_fastqc.zip"
    log:
        "output/logs/fastqc/{prefix}.log"
    shell:
        "fastqc {input} -o output/fastqc > {log}"

#bwa mem threads
#rule bwa_mem_threads:
#    input:
#        genome="data/bwa_index/sacCer3.fa",
#        fq1="data/samples/{sra}_1.fastq",
#        fq2="data/samples/{sra}_2.fastq"
#    output:
#        "output/mapped/{sra}.sam"
#    threads: 
#        4
#    log:
#        "output/logs/bwa/{sra}.log"
#    shell:
#        "bwa mem -t {threads} {input.genome} {input.fq1} {input.fq2} 1> {output} 2> {log}"

#bwa mem param
rule bwa_mem_param:
    input:
        genome="data/bwa_index/sacCer3.fa",
        fq1="data/samples/{sra}_1.fastq",
        fq2="data/samples/{sra}_2.fastq"
    output:
        "output/mapped/{sra}.sam"
    threads: 
        4
    params:
        extra="-v 1"
    log:
        "output/logs/bwa/{sra}.log"
    shell:
        "bwa mem -t {threads} {params.extra} {input.genome} {input.fq1} {input.fq2} 1> {output} 2> {log}"

#practice: sam to bam
#rule sam_to_bam:
#    input:
#    output:
#    shell:

#samtools merge
sras = ["SRR922442", "SRR922443"]
rule merge_bam:
    input:
        expand("output/mapped/{sra}.bam", sra=sras)
    output:
        "output/mapped/final.merge.bam"
    threads:
        4
    params:
        extra=""
    shell:
        "samtools merge -@ {threads} {params.extra} {output} {input}"

# python script integration
rule script:
    input:
        "data/input.tsv"
    output:
        "output/figures/script.png"
    script:
        "scripts/example.py"

# jupyter notebook
rule jupyter:
    input:
        "data/input.tsv"
    output:
        "output/figures/jupyter.png"
    log:
        # optional path to the processed notebook
        notebook="logs/notebooks/processed_notebook.ipynb"
    notebook:
        "notebooks/test.ipynb"