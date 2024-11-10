#!/bin/bash

#Author: ZHOU YUQI

#This Script facilitate the mapping of sequences to the reference genome
#We assume that the sequences' quality have been assessed and trimmed
#This script contains the use of bwa and samtools. 
#Therefore, to successfully run the scripts, bwa and samtools are required to be installed as binary files at /bin folder under the conda env or your local folder.

#The script takes reference genome, trimmed fastq files as input arguments
#The input order should be: referece genome, trimmed fastq file 1, trimmed fastq file 2 (Paired-end reads mapping)
#Or reference genome, trimmed fastq file (single end reads mapping)
#Script can solve the single/paired situation automatically.

#Input format: the script is designed for dealing with the output files of run_trimmomatic.sh
#Therefore, the input file name format is : identifier.fasta, identifier_1.trim.fastq, identifier_2.trim.fastq

#Absolute path of input files are required to ensure successful running.
#output files will be in a folded named 'results' at current location

#help function to remind the input format
show_help() {
  echo "
  Format to run run_mapping.sh:
    ./run_mapping.sh <reference genome> <trimmed file 1> <trimmed file 2>
    or:
    ./run_mapping.sh <reference genome> <trimmed file> 

  Required input filename format:
    <reference genome>: identifier.fasta
    <trimmed file 1>: identifier_1.trim.fastq
    <trimmed file 2>: identifier_2.trim.fastq
    <trimmed file>: identifier.trim.fastq

  For more details, please view the script comments.
     
  "
}

#Condition to give the help infor
if [ "$1" = -h ] || [ "$1" = --help ]
then 
    show_help
    exit 0
fi

SRA=$1
cd ./$SRA/ref_gene
ref_genome=($(ls *fasta)) #find ref_genome

#index the reference genome
bwa index $ref_genome

#find trimmed gz files
cd ../raw_data/trimmed_fastq
mapping_reads_1_gz=($(ls *_1.trim.fastq.gz))
mapping_reads_2_gz=($(ls *_2.trim.fastq.gz))
gunzip ./$mapping_reads_1_gz
gunzip ./$mapping_reads_2_gz

#find trimmed ungzip files
mapping_reads_1=($(ls *_1.trim.fastq))
mapping_reads_2=($(ls *_2.trim.fastq))



cd ..
cd ..

mkdir -p ./mapping/sam

#Now position is under SRA
out1=${SRA}.aligned.sam
bwa mem ./ref_gen/$ref_genome ./raw_data/trimmed_fastq/$mapping_reads_1 ./raw_data/trimmed_fastq/$mapping_reads_2 > ./mapping/sam/$out1

#Convert SAM to BAM with samtools
mkdir -p ./mapping/bam
out2=${SRA}.aligned.bam
samtools view -S -b ./mapping/sam/$out1 > ./mapping/bam/$out2

#Sort the BAM file with samtools
out3=${SRA}.aligned.sorted.bam
samtools sort -o ./mapping/bam/$out3 ./mapping/bam/$out2

#Index the sorted BAM file
samtools index ./mapping/bam/$out3

echo "The mapping is done successfully."






