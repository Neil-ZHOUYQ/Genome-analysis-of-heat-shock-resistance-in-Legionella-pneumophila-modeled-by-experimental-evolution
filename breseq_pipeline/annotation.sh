#!/bin/bash
# Author: Jason
# Date: 2024/12/11
# Set input file and target directory

sorted_file="sorted_SRArun.txt"
vcf_folder="/hdd1/home/f24_wyli/project/pipeline/vc_fastp/grouped_vcf"  # Folder containing grouped VCF files
output_folder="annotated_vcf"  # Folder to store annotated VCF files
log_file="snpeff_annotation.log"  # Log file

# Set SnpEff parameters
snpeff_jar="/hdd1/home/f24_wyli/project/pipeline/tools/snpEff/snpEff.jar"  # Path to SnpEff
genome="L_p"       # Genome database to use

# Create output directory and log file
mkdir -p "$output_folder"
echo "SnpEff Annotation Log - $(date)" > "$log_file"

# Loop through the grouped folders and find each group's VCF files
for group_folder in "$vcf_folder"/*; do
    group_name=$(basename "$group_folder")
    annotated_group_folder="$output_folder/$group_name"
    mkdir -p "$annotated_group_folder"

    echo "Processing group: $group_name" | tee -a "$log_file"

    for vcf_file in "$group_folder"/*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            vcf_basename=$(basename "$vcf_file")
            annotated_vcf="$annotated_group_folder/$vcf_basename"

            echo "Annotating file: $vcf_file" | tee -a "$log_file"
            echo "Running command: java -Xmx4g -jar $snpeff_jar -v $genome $vcf_file > $annotated_vcf" | tee -a "$log_file"
            
            # Run SnpEff and log output
            java -Xmx4g -jar "$snpeff_jar" "$genome" "$vcf_file" > "$annotated_vcf" 2>>"$log_file"

            if [[ $? -eq 0 ]]; then
                echo "Annotation completed: $annotated_vcf" | tee -a "$log_file"
            else
                echo "Annotation failed: $vcf_file" | tee -a "$log_file"
            fi
        else
            echo "Warning: VCF file not found $vcf_file" | tee -a "$log_file"
        fi
    done
done

echo "All files annotated - $(date)" | tee -a "$log_file"
