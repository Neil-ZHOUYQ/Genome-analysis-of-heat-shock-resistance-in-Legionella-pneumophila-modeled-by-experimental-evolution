SRA=$1

mkdir -p ./$SRA

cd ./$SRA

mkdir -p ./ref_gen/ ./mapping ./raw_data ./calling

cd ..



cp ./reference/* ./$SRA/ref_gen/



./codes/download.sh $SRA
./codes/fastqc_trim.sh $SRA
./codes/run_mapping.sh $SRA
./codes/run_variant_calling.sh $SRA
./codes/annotation.sh $SRA
