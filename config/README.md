Here about how to configure your data

The config.yaml is a file where you must specify the following information:
  - genome: here you must write the path where you store the genome fasta file.

    For example
      genome: "/home/user/data/genome.fa"
  
  - samples: here you must write the root of your pairends.

      For example: If the samples are Sample01_1.fasta and Sample01_2.fasta, the root will be: Sample01
        samples: Sample01

      If you want to study more than one sample, the root must be inside brackets.
        samples: ["Sample01","Sample02","Sample03"]
  
  - units: This is the number of each pairends samples.
           units: ["1","2"]

  - separation: this is the symbol that is between the root and the unit in the sample name. 
           separation: _

However, the samples you are going to study, must be inside the data/sample directory from this github repository. 

    
