# barley_variant_calling
Repository practice project in EEAD-CSIC

## INTRODUCTION

The main objective is to analyse the variant calling with GATK. To achive that goal I will follow the following process.
1. Prepare Data
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

To configure this workflow to use your own data, go to config directory and follow the instructions. 

# TO USE THIS REPOSITORY
To use all the program at one, write in the TErminal the following code:


    snakemake -p -c4

-c or --cores you can select the number of cores you want

# STEP BY STEP
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
#### 3.2.1 Map the paired ends reads and get the output BAM file
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

### 3.4. Visualize the map reads in IGV

#Step 5: Depth calculation

rule depth_calc:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/depth/{sample}_depth.csv"
    shell:
        "samtools depth {input} > {output}"

#Step 6: Show th depth in a plot

rule plot_depth:
    input:
        "results/depth/{sample}_depth.csv"
    output:
        "results/depth/plots/{sample}.svg"
    script:
        "plot-depth.py"

        
## 4. FIGURE MAPPING CHROMOSOME
### 4.1. Depth calculation

First of all, we need to calculate the depth. For this calculation, we use the rule depth_calc, where we get a .csv output. 

To run the rule, write this code in the terminal:

    snakemake -p results/depth/A_1_20_depth.csv  -c4

### 4.1. Show depth in a plot
For this step, we are going to use python 3. However, you can use the programm you want to visualize the data.

You can find the program I did in the directory workflow/scripts

## 5. VARIANT CALLING

Software used in this practical:
- Picard 
- GATK
  
Installing GATK


### 5.1. Mark Duplicates
We are using the 

### 5.2. Haplotype Calling![Uploading dag-plot-included.svg…]()

GATK function



## TO VISUALIZE THE PROCESS WITH DAG

![dag-plot-included](https://github.com/carmenmiravete/barley_variant_calling/assets/151924636/67e738ff-4309-4036-9ccf-b32379afdba8)


    snakemake --dag -p | dot -Tsvg > dag.svg^C

 
