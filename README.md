# G5pipe & Genomic variation analysis of Legionella pneumophila with heat-resistance 
The repository keeps the codes of a bioinformatic research project, tagetting to analysis genomic variations  in Legionella pneumophila with heat resistance modeled by experimental evolution, in The Chinese University of Hong Kong. In the project, an analysis tool G5pipe was developed, displaying significant capability of  finding SNPs. The project focused on anlyzing the SNPs found by both G5pipe and breseq(http://barricklab.org/breseq).
## G5pipe_structure

![G5pipe_structure](images/G5pipe_structure.png)

G5pipe has scripts for each step, which facilitates customizing workflow and debugging. 

## Preparation for run G5pipe
1. A conda environment containing trimmomatic and snpEff
2. A snpEff database constructed by yourself. A reference blog link is: https://www.cnblogs.com/JewelZ/p/17364098.html.
3. A genome.fasta and a genome.gff

## Run G5pipe without annotation
snpEff, snpEff database and genome.gff are not needed when annotation is not required.
2 name-specified directories are needed under current directory: codes and reference, which contains all shell scripts files and genome FASTA file respectively. Results are structured in different directories 
[picture]

## Run G5pipe with annotation
The species name in annotation.sh should be modified to the species where sequencing data comes from. My example used Legionella pneumophila as reference genome.
2 name-specified directories are needed under current directory: codes and reference. ./codes contains all shell scripts files. ./reference contains genome.fasta and genome.gff. Results are structured in different directories 
[picture]

## Run G5pipe using our project dataset
[picture]


