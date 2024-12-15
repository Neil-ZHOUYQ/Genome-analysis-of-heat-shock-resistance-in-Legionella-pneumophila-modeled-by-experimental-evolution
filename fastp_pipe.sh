#!/bin/bash
# This script uses fastp tool to do the data clean job
# Author: Jason
# Date: 2024/12/11

# Directories Args
input_dir="raw_data"  # Input raw data folder
output_dir="output_fastq_files"  # Output FASTQ folder
html_report_dir="qc_reports"  # Folder for HTML reports
json_report_dir="json_reports"  # Folder for JSON reports

# Create output directories
mkdir -p "$output_dir" "$html_report_dir" "$json_report_dir"

# Loop through all paired SRR subfolders (ID range from SRR24205963 to SRR24206065)
for i in {24205963..24206065}; do
    # Construct SRR folder path
    sr_folder="$input_dir/SRR${i}"

    # Check if SRR subfolder exists
    if [ ! -d "$sr_folder" ]; then
        echo "Warning: $sr_folder does not exist. Skipping..."
        continue
    fi

    # Check if paired input files exist (_1.fastq and _2.fastq)
    input_r1="${sr_folder}/SRR${i}_1.fastq"
    input_r2="${sr_folder}/SRR${i}_2.fastq"

    if [ ! -f "$input_r1" ] || [ ! -f "$input_r2" ]; then
        echo "Warning: One or both input files for SRR${i} are missing. Skipping..."
        continue
    fi

    # Construct output file paths
    output_r1="$output_dir/SRR${i}_1_clean.fastq"
    output_r2="$output_dir/SRR${i}_2_clean.fastq"
    html_report="$html_report_dir/SRR${i}_qc_report.html"
    json_report="$json_report_dir/SRR${i}_qc_report.json"
    log_file="$output_dir/SRR${i}_fastp.log"  # Log file

    # Run fastp for quality control
    fastp -i "$input_r1" -I "$input_r2" -o "$output_r1" -O "$output_r2" \
          --html "$html_report" --json "$json_report" \
          --cut_mean_quality 30 --length_required 50 \
          --adapter_sequence AGATCGGAAGAGC --adapter_sequence_r2 AGATCGGAAGAGC \
          -w 4 &> "$log_file"  # Use 4 threads and redirect output to log file

    # Output log
    if [ $? -eq 0 ]; then
        echo "Processing sample SRR${i} completed successfully."
    else
        echo "Error processing sample SRR${i}. Check log file $log_file for details."
    fi
done

echo "Quality control for all samples completed!"
