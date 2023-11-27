rule map_reads:
    input:
        "data/GCA_904849725_genome.fa,
        "data/samples/{sample}.fastq"
    output:
        "results/mapped/{sample}.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa mem {input} | samtools view -b - > {output}"

