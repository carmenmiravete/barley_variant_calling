# barley_variant_calling
Repository practice project in EEAD-CSIC

## Introduction

The main objective is to analyse the variant calling with GATK. To achive that goal I will follow the following process.

1. Quality Control
2. Preprocessing
3. (Assembly)
4. Mapping
5. (Annotation)
6. Variant calling


Genome reference MorexV
Samples: 32 samples of barley

## STARTING WITH SNAKEMAKE
This tutorial introduces the text-based workflow system Snakemake: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html 

#Alll the rules are in the Snakemakefile in workflow Directory
## 1. DESCOMPRESS FILES
It is possible that the files are compressed.

Make a rule to decompress the files. 

In the terminal: 

´snakemake -p data/samples/{A_1_21_1,A_1_21_2}.fastq -c4´

    snakemake -p data/samples/{A_1_21_1,A_1_21_2}.fastq -c4



## 2. QUALITY CONTROL

I will start the project doing a QUality Control with FastQC. These sequences must be verified for quality control, to ensure that the raw data does not have problems that could affect the biological analysis. 

Process:
  Import of data from FastQ files 2 by 2. 
  The program will provide a overview of the data
  Summary graphs
  Export results to an HTML based permanent report

  For this process you must use 
