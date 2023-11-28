# **PROJECT CSIC PRACTICE**

#Introduction

The main objective is to analyse the variant calling with GATK. To achive that goal I will follow the following process.

1. Quality Control
2. Preprocessing
3. (Assembly)
4. Mapping
5. (Annotation)
6. Variant calling


Genome reference MorexV
Samples: 32 samples of barley

## 0. DESCOMPRESS FILES
It is possible that the files are compressed.

I made a rule descompressed:

In the terminal: snakemake -p data/samples/{A_1_21_1,A_1_21_2}.fastq -c4

## 1. QUALITY CONTROL

I will start the project doing a QUality Control with FastQC. These sequences must be verified for quality control, to ensure that the raw data does not have problems that could affect the biological analysis. 

Process:
  Import of data from FastQ files 2 by 2. 
  The program will provide a overview of the data
  Summary graphs
  Export results to an HTML based permanent report
  
