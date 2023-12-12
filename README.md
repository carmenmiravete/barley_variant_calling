# barley_variant_calling
Repository practice project in EEAD-CSIC

## Introduction

The main objective is to analyse the variant calling with GATK. To achive that goal I will follow the following process.
1. Prepare DAta
1. Quality Control
2. Mapping
3. Variant calling


Genome reference MorexV
Samples: 32 samples of barley

## STARTING WITH SNAKEMAKE
This tutorial introduces the text-based workflow system Snakemake: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html 

Alll the rules are in the Snakemakefile in workflow Directory

## DATA

We need a paired end couple of samples or a single end. In this repository, we have two samples "A_1_20_1" and "A_1_20_2" in the data/samples directory. Besides, in the data directory we can find the whole genome reference: GCA_904849725_genome.fa


## 1. DESCOMPRESS FILES
It is possible that the files are compressed.

Make a rule to decompress the files. 

In the terminal: 

    snakemake -p data/samples/{A_1_20_1,A_1_20_2}.fastq.gz -c4


## 2. QUALITY CONTROL

I will start the project doing a QUality Control with FastQC. These sequences must be verified for quality control, to ensure that the raw data does not have problems that could affect the biological analysis. 

Process:
  Import of data from FastQ files 2 by 2. 
  The program will provide a overview of the data
  Summary graphs
  Export results to an HTML based permanent report

 
    snakemake -p results/reports/data/samples/{A_1_20_1,A_1_20_1}.html  -c4

## 3. MAPPING
### Files in use:
- GCA_904849725_genome.fa
- SAMPLES: A_1_20_1 and A_1_20_2 

### 3.1. Index the reference genome 

It can be done from Terminal:

    bwa index -a bwtsw data/GCA_904849725_genome.fa

Or you can find the respective rule in the Snakefile called rule ref_genome.
### 3.2. Mapping (bwa mem)
### 3.2.1 Map the paired ends reads and get the output BAM file
The mapping reads information is in the rule map_reads in the Snakefile. 

To run this rule, introduce the following code in the terminal

    snakemake -p results/sorted/A_1_21_RG.bam

#### 3.2.1. Map the single end reads 
It can be the case of map single end reads.
Single end Read alignment:
    bwa aln -f data/GCA_904849725_genome.fa Single_End_Sample.fastq 

### 3.3. Sort the bam file

There ir a rule called rule sort_bam. To run this rule, we must write the following code in the Terminal:

    snakemake -p results/sorted/A_1_20_sorted.bam  -c4

## 






  
