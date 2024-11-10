#!/bin/bash

#Author: ZHOU YUQI

#This script facilitates the variant calling analysis of run_mapping.sh. 
#This script is designed for dealing with the output files of run_mapping.sh
#The script contains the use of bcftools.
#Therefor, to successfully run the scripts, btftools is required to be installed 
#as binary files at /bin folder under the conda env or your local folder.

#The script takes reference genome, aligned bam files as input arguments.
#The input order should be: referece genome, aligned bam files.
#The input format is: ./run_variant.sh <reference genome> <aligned file> 

#Absolute path of input files are required to ensure successful running.
#output files will be in a folded named 'results' at current location



SRA=$1

#find aligned_file and ref_genome
cd ./$SRA/mapping/bam
aligned_file=($(ls *aligned.sorted.bam))

cd ..
cd ..
cd ../ref_gen
ref_genome=($(ls *fasta))
cd ..


#Mark duplicated reads
out1=${SRA}_raw.bcf
mkdir -p calling/bcf
bcftools mpileup -O b -o calling/bcf/$out1 -f ./ref_gene/$ref_genome ./mapping/bam/$aligned_file

#Identify SNPs using bcftools
out2=${SRA}_variants.vcf
mkdir -p calling/vcf
bcftools call --ploidy 1 -m -v -o calling/vcf/$out2 calling/bcf/$out1

#Filter SNPs
out3=${SRA}_final_variants.vcf
vcfutils.pl varFilter calling/vcf/$out2 > calling/vcf/$out3 

echo "variant calling is done successfully."



