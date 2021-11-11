# first example
rule bwa_map:
    input:
        "data/chr1.fa",
        "data/samples/A.fastq"
    output:
        "output/mapped/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

# example of wildcard (generalized rule)
rule bwa_map_generalize:
    input:
        "data/chr1.fa",
        "data/samples/{sample}.fastq"
    output:
        "output/mapped/{sample}.wildcards.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

# threads
rule bwa_map_threads:
    input:
        "data/chr1.fa",
        "data/samples/{sample}.fastq"
    output:
        "output/mapped/{sample}.threads.bam"
    threads: 4
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

# params:
rule bwa_map_params:
    input:
        "data/chr1.fa",
        "data/samples/{sample}.fastq"
    output:
        "output/mapped/{sample}.params.bam"
    threads: 4
    params: "-v 1"
    shell:
        "bwa mem -t {threads} {params} {input} | samtools view -Sb - > {output}"

# example of expand function
SAMPLES = ["A", "B"]
rule samtools_merge_expand:
    input:
        expand("output/mapped/{sample}.params.bam", sample = SAMPLES)
    output:
        "output/merge/merge.expand.bam"
    shell:
        "samtools merge {output} {input}"

# example of glob_wildcards function
rule samtools_merge_glob_wildcards:
    input:
        expand("output/mapped/{sample}.params.bam", sample = glob_wildcards("data/samples/{sample, [AB]}.fastq").sample)
    output:
        "output/merge/merge.glob_wildcards.bam"
    shell:
        "samtools merge {output} {input}"

# conda environment
rule samtools_merge_conda:
    input:
        expand("output/mapped/{sample}.params.bam", sample = SAMPLES)
    output:
        "output/merge/merge.conda.bam"
    conda:
        "environment.yaml"
    shell:
        "samtools merge {output} {input}"

# customized function using wildcards
import os
def find_input(wildcards):
    fs = []
    for f in os.listdir("data/samples/"):
        if "fastq" in f:
            fs.append("output/mapped/" + f.replace(".fastq", ".params.bam"))
    return fs

rule samtools_merge_customized:
    input:
        find_input
    output:
        "output/merge/merge.bam"
    shell:
        "samtools merge {output} {input}"
    
