#!/bin/bash

# GNBF5030 Project Analysis Pipeline (breseq-leading version)
# This script performs data analysis including download, QC, preprocessing, and variant calling
# Author: Jason
# Date: 2024/11/22

# Exit on error
set -e

# Enable command printing for debugging
set -x

# Error handling function with more detailed output
error_handler() {
    echo "Error occurred in script at line: $1"
    echo "Command that failed: $BASH_COMMAND"
    echo "Working directory: $(pwd)"
    exit 1
}

trap 'error_handler ${LINENO}' ERR

# Define base directory and create directory structure
BASE_DIR="$(pwd)/analysis_results"
mkdir -p "${BASE_DIR}"/{raw_data,reference,results,logs}
cd "${BASE_DIR}"

# Checkpoint file to track progress
CHECKPOINT_FILE="${BASE_DIR}/logs/checkpoint.log"
PROGRESS_FILE="${BASE_DIR}/logs/progress.log"

# Create temporary SRR list file
cat > srr_list.txt << 'EOL'
SRR24205963
SRR24205964
SRR24205965
SRR24205966
SRR24205967
SRR24205968
SRR24205969
SRR24205970
SRR24205971
SRR24205972
SRR24205973
SRR24205974
SRR24205975
SRR24205976
SRR24205977
SRR24205978
SRR24205979
SRR24205980
SRR24205981
SRR24205982
SRR24205983
SRR24205984
SRR24205985
SRR24205986
SRR24205987
SRR24205988
SRR24205989
SRR24205990
SRR24205991
SRR24205992
SRR24205993
SRR24205994
SRR24205995
SRR24205996
SRR24205997
SRR24205998
SRR24205999
SRR24206000
SRR24206001
SRR24206002
SRR24206003
SRR24206004
SRR24206005
SRR24206006
SRR24206007
SRR24206008
SRR24206009
SRR24206010
SRR24206011
SRR24206012
SRR24206013
SRR24206014
SRR24206015
SRR24206016
SRR24206017
SRR24206018
SRR24206019
SRR24206020
SRR24206021
SRR24206022
SRR24206023
SRR24206024
SRR24206025
SRR24206026
SRR24206027
SRR24206028
SRR24206029
SRR24206030
SRR24206031
SRR24206032
SRR24206033
SRR24206034
SRR24206035
SRR24206036
SRR24206037
SRR24206038
SRR24206039
SRR24206040
SRR24206041
SRR24206042
SRR24206043
SRR24206044
SRR24206045
SRR24206046
SRR24206047
SRR24206048
SRR24206049
SRR24206050
SRR24206051
SRR24206052
SRR24206053
SRR24206054
SRR24206055
SRR24206056
SRR24206057
SRR24206058
SRR24206059
SRR24206060
SRR24206061
SRR24206062
SRR24206063
SRR24206064
SRR24206065
EOL

# Initialize checkpoint file if it doesn't exist
if [ ! -f "${CHECKPOINT_FILE}" ]; then
    touch "${CHECKPOINT_FILE}"
fi

# Function to check if a sample has been completed
is_sample_completed() {
    local srr=$1
    grep -q "^${srr}:COMPLETED$" "${CHECKPOINT_FILE}"
    return $?
}

# Function to mark a sample as completed
mark_sample_completed() {
    local srr=$1
    echo "${srr}:COMPLETED" >> "${CHECKPOINT_FILE}"
}

# Function to log progress
log_progress() {
    local message=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${message}" >> "${PROGRESS_FILE}"
    echo "${message}"
}

# Download reference genome and annotation files
log_progress "Downloading reference files..."
cd reference
if [ ! -f "GCF_001752745.1_ASM175274v1_genomic.gbff" ]; then
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/752/745/GCF_001752745.1_ASM175274v1/GCF_001752745.1_ASM175274v1_genomic.gbff.gz
    gunzip *.gz
fi
cd ..

# Check if reference file exists and is readable
if [ ! -r "reference/GCF_001752745.1_ASM175274v1_genomic.gbff" ]; then
    log_progress "Error: Reference file not found or not readable"
    exit 1
fi

# Process each SRR file sequentially
log_progress "Starting sequential processing of samples..."
while read -r srr; do
    # Skip if already completed
    if is_sample_completed "${srr}"; then
        log_progress "Skipping ${srr} (already completed)"
        continue
    fi

    log_progress "Processing ${srr}"
    
    # Create directories for current SRR
    raw_dir="${BASE_DIR}/raw_data/${srr}"
    results_dir="${BASE_DIR}/results/${srr}"
    mkdir -p "${raw_dir}" "${results_dir}"
    
    # Download SRA file (if not already downloaded)
    if [ ! -d "${raw_dir}" ] || [ -z "$(ls -A ${raw_dir})" ]; then
        log_progress "Downloading ${srr}..."
        prefetch "${srr}"
    else
        log_progress "SRA file for ${srr} already exists, skipping download"
    fi
    
    # Convert to FASTQ format (if not already converted)
    if ! ls "${raw_dir}"/*.fastq >/dev/null 2>&1; then
        log_progress "Converting ${srr} to FASTQ..."
        if ! fastq-dump --split-files --outdir "${raw_dir}" "${srr}" 2>"${raw_dir}/fastq-dump.log"; then
            log_progress "Error in fastq-dump for ${srr}"
            cat "${raw_dir}/fastq-dump.log"
            continue
        fi
    else
        log_progress "FASTQ files for ${srr} already exist, skipping conversion"
    fi
    
    # Check if FASTQ files were created
    if ! ls "${raw_dir}"/*.fastq >/dev/null 2>&1; then
        log_progress "Error: No FASTQ files found for ${srr}"
        continue
    fi
    
    # Get the number of available CPU cores and set the number of threads dynamically
    THREADS=$(nproc)

    # Defines the path to the flag file
    FLAG_FILE="${results_dir}/breseq_completed.flag"

    # Check if the input file exists to avoid long argument lists
    INPUT_FILES=(${raw_dir}/*.fastq)
    if [ ${#INPUT_FILES[@]} -eq 0 ]; then
        log_progress "No FASTQ files found in ${raw_dir} for ${srr}, skipping analysis."
        exit 1
    fi

    # run breseq（if have not completed）
    if [ ! -f "${FLAG_FILE}" ]; then
        log_progress "Running breseq for ${srr}..."

        if ! breseq \
            -r "${BASE_DIR}/reference/GCF_001752745.1_ASM175274v1_genomic.gbff" "${INPUT_FILES[@]}" \
            -o "${results_dir}" \
            -j "${THREADS}" \
            &> "${results_dir}/breseq.log"; then
            log_progress "Error in breseq for ${srr}"
            cat "${results_dir}/breseq.log"

            # Clean up possible incomplete output
            rm -rf "${results_dir}"

            # If in a loop, use continue. Otherwise, use exit 1
            # continue  
            exit 1      
        else
            # Create the flag file On successful completion
            touch "${FLAG_FILE}"
        fi
    else
        log_progress "Breseq results for ${srr} already exist, skipping analysis."
    fi

    # Optional: Remove FASTQ files to save space after successful processing
    # if [ -f "${results_dir}/output/index.html" ]; then
    #     rm "${raw_dir}"/*.fastq
    # fi
    
    mark_sample_completed "${srr}"
    log_progress "Completed processing ${srr}"
    
done < srr_list.txt

# Generate final report
log_progress "Generating final report..."
total_samples=0
successful_samples=0
failed_samples=0

while read -r srr; do
    total_samples=$((total_samples + 1))
    if [ -f "${BASE_DIR}/results/${srr}/output/index.html" ]; then
        successful_samples=$((successful_samples + 1))
        log_progress "${srr}: Success"
    else
        failed_samples=$((failed_samples + 1))
        log_progress "${srr}: Failed"
        if [ -f "${BASE_DIR}/results/${srr}/breseq.log" ]; then
            cat "${BASE_DIR}/results/${srr}/breseq.log" >> "${PROGRESS_FILE}"
        fi
    fi
done < srr_list.txt

# Final status report
log_progress "Pipeline completed!"
log_progress "Total samples: ${total_samples}"
log_progress "Successful: ${successful_samples}"
log_progress "Failed: ${failed_samples}"
log_progress "Detailed logs available in: ${PROGRESS_FILE}"

# Clean the temporary file 
# rm ${srr}
# Print system information
echo "System information:"
df -h . >> "${PROGRESS_FILE}"

# Optional cleanup
# rm srr_list.txt
