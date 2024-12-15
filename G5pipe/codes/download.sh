SRA=$1

# mkdir -p ./$SRA

cd ./$SRA

# mkdir -p ./ref_gen/ ./mapping ./raw_data ./calling

cd ./raw_data

prefetch $SRA

cp ./$SRA/${SRA}.sra .
rm -r ./$SRA

echo 'The .sra file will be downloaded is:'
ls

echo '$SRA is downloading...'
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ./${SRA}.sra

echo '$SRA downloaded files are:'
ls ./fastq


cd ..
cd ..