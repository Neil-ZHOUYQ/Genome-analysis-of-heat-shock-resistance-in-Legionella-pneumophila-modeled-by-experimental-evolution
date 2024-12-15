#!/bin/bash
# Author: Jason
# Date: 2024/12/11
# Input folder and reference genome path
input_dir="/hdd1/home/f24_wyli/project/pipeline/fastp_quality_control/output_fastq_files"  # Folder storing cleaned FASTQ files
reference_genome="/hdd1/home/f24_wyli/project/pipeline/analysis_results/reference/GCF_001752745.1_ASM175274v1_genomic.gbff"  # Path to reference genome
output_base_dir="breseq_outputs"  # Base folder to save breseq output results

# Create the output base folder (if it doesn't exist)
mkdir -p "$output_base_dir"

# Loop through all SRR subfolders
for folder in "$input_dir"/*/; do
    # Get SRR ID (folder name)
    srr_id=$(basename "$folder")

    # Input files (_1_clean.fastq and _2_clean.fastq)
    input_r1="$folder/${srr_id}_1_clean.fastq"
    input_r2="$folder/${srr_id}_2_clean.fastq"

    # Create breseq output folder for the corresponding SRR ID
    output_dir="$output_base_dir/$srr_id"
    mkdir -p "$output_dir"

    # Run breseq
    breseq -r "$reference_genome" "$input_r1" "$input_r2" -o "$output_dir"

    echo "Variant calling for $srr_id completed. Results saved in $output_dir."
done

echo "Breseq variant calling for all samples completed!"
