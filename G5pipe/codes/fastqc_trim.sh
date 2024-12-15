SRA=$1
cd ./$SRA/raw_data


mkdir fastqc_results
fastqc ./fastq/*.gz -o ./fastqc_results

cd fastqc_results
for file in *.zip; do unzip $file; done

cd ..
mkdir ./trimmed_fastq
cd fastq 

files=($(ls))
in1=${files[0]}
in2=${files[1]}

out1=${SRA}_1.trim.fastq.gz
out2=${SRA}_1un.trim.fastq.gz
out3=${SRA}_2.trim.fastq.gz
out4=${SRA}_2un.trim.fastq.gz
trimmomatic PE $in1 $in2 ../trimmed_fastq/$out1 ../trimmed_fastq/$out2 ../trimmed_fastq/$out3 ../trimmed_fastq/$out4 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:25

cd ..
cd ..
cd ..