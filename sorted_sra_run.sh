#!/bin/bash

# Set input file and target directory
sorted_file="sorted_SRArun.txt"
vcf_folder="/hdd1/home/f24_wyli/project/pipeline/vc_fastp/breseq_outputs"  # Folder storing all VCF files
output_folder="grouped_vcf"

# Create output folder
mkdir -p "$output_folder"

# Define the variant location to remove
chromosome="NZ_CP015927"
position="173322"

# Iterate through sorted_SRArun.txt and create folders by group
while read -r line; do
    # Extract SRR ID and group name
    srr_id=$(echo "$line" | awk '{print $1}')
    group_name=$(echo "$line" | awk '{print $2}')
    
    # Create corresponding group folder
    group_folder="$output_folder/$group_name"
    mkdir -p "$group_folder"
    
    # Process the corresponding VCF file
    vcf_file="$vcf_folder/${srr_id}.vcf"
    if [[ -f "$vcf_file" ]]; then
        # Remove specific variant and save to the target folder
        grep -v -P "^${chromosome}\t${position}\t" "$vcf_file" > "$group_folder/${srr_id}.vcf"
    else
        echo "Warning: VCF file $vcf_file not found"
    fi
done < "$sorted_file"

echo "File processing and grouping completed!"
