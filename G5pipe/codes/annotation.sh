SRA=$1
cd ./${SRA}/calling/vcf

java -jar ~/apps/snpEff/snpEff.jar -c ~/apps/snpEff/snpEff.config Philadelphia_1_ATCC ./${SRA}_final_variants.vcf > ./annotated_${SRA}.vcf

echo "VCF file is annotated successfully."

cd ..
cd ..
cd ..
