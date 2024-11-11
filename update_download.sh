#!/bin/bash

################################################################################
# Script: download.sh 
# Description: Download and process SRA data for L. pneumophila sequencing analysis
#
# Usage: ./download.sh <SRA_ID> <REF_GENOME_SITE>
#
# Arguments:
#   SRA_ID          - The SRA accession number to download
#   REF_GENOME_SITE - Reference genome location
#
# Options:
#   -h, --help    - Display this help message
#
# Output:
#   - Downloaded SRA files in ./${SRA}/raw_data/
#   - Converted FASTQ files in ./${SRA}/raw_data/fastq/
#   - Log files in ./${SRA}/logs/
#
# Directory Structure:
#   ./${SRA}/
#   ├── raw_data/   - Contains raw sequencing data
#   ├── ref_gen/    - Reference genome files
#   ├── mapping/    - Mapping results
#   ├── calling/    - Variant calling results
#   └── logs/       - Log files
#
# Dependencies:
#   - SRA Toolkit (prefetch, fastq-dump)
#   - md5sum
#
# Example:
#   ./download.sh SRR12345678 /path/to/reference
#
# Author: Neil & Jason
# Date: 11/11/2024
# Version: 1.0
#
# Notes:
#   - Requires sufficient disk space for SRA and FASTQ files
#   - Internet connection required for SRA download
#   - Paired-end sequencing data is assumed
################################################################################

# Exit on any error
set -e

# Function to display help message
show_help() {
    head -n 32 "$0" | grep "^#" | sed 's/^#//'
    exit 0
}

# Parse command line arguments for help
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    show_help
fi

# Set up logging function
setup_logging() {
    local log_dir="./${SRA}/logs"
    mkdir -p ${log_dir}
    # Redirect stdout and stderr to both console and log files
    exec 1> >(tee -a "${log_dir}/${SRA}_download.log")
    exec 2> >(tee -a "${log_dir}/${SRA}_download.error.log" >&2)
    echo "Started download process at $(date)"
}

# Set the path
# SRA_PATH="./${SRA}/raw_data/${SRA}/${SRA}.sra"

# Function to check required tools
check_dependencies() {
    local tools=("prefetch" "fastq-dump" "md5sum")
    for tool in "${tools[@]}"; do
        if ! command -v $tool &> /dev/null; then
            echo "Error: Required tool '$tool' is not installed."
            echo "Please install SRA Toolkit for prefetch and fastq-dump"
            echo "Installation options:"
            echo "  conda install -c bioconda sra-tools"
            echo "  or visit: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit"
            exit 1
        fi
    done
    echo "All required tools are available"
}

# Function to validate input parameters
validate_inputs() {
    if [ $# -ne 2 ]; then
        echo "Error: Invalid number of parameters"
        echo "Usage: $0 <SRA_ID> <REF_GENOME_SITE>"
        echo "Try '$0 --help' for more information"
        exit 1
    fi
    
    # Validate SRA ID format
    if ! [[ $1 =~ ^[SDE]RR[0-9]{6,} ]]; then
        echo "Warning: SRA ID format may be incorrect. Expected format: SRR######"
    fi
}

# Function to create directory structure
create_directories() {
    echo "Creating directory structure..."
    for dir in ref_gen mapping raw_data calling logs; do
        mkdir -p ./${SRA}/${dir}
        echo "Created directory: ${SRA}/${dir}"
    done
}

# Function to download SRA data
download_sra() {
    echo "Downloading SRA file: ${SRA}"
    if prefetch ${SRA}; then
        echo "SRA download completed successfully"
    else
        echo "Error: Failed to download SRA file"
        exit 1
    fi
}

# Function to convert SRA to FASTQ format
convert_to_fastq() {
    echo "Converting SRA to FASTQ format..."
    mkdir -p fastq
    
    # Run fastq-dump with appropriate parameters
    fastq-dump --outdir fastq \
               --gzip \
               --skip-technical \
               --readids \
               --dumpbase \
               --split-3 \
               --clip \
               ${SRA}/${SRA}.sra
    
    if [ $? -eq 0 ]; then
        echo "FASTQ conversion completed successfully"
    else
        echo "Error: FASTQ conversion failed"
        exit 1
    fi
}

# Function to verify output files
verify_files() {
    echo "Verifying file integrity..."
    local fastq_dir="./fastq"
    
    # Check for paired-end files
    for end in 1 2; do
        if [ ! -f "${fastq_dir}/${SRA}_${end}.fastq.gz" ]; then
            echo "Error: Missing FASTQ file: ${SRA}_${end}.fastq.gz"
            exit 1
        fi
    done
    
    # Generate checksums
    md5sum ${fastq_dir}/${SRA}_*.fastq.gz > ${fastq_dir}/checksums.md5
    echo "File verification completed successfully"
}

# Function to clean up temporary files
cleanup() {
    echo "Cleaning up temporary files..."
    if [ -f "${SRA}/${SRA}.sra" ]; then
        rm -f ${SRA}/${SRA}.sra
        echo "Removed temporary SRA file"
    fi
}

# Main execution function
main() {
    # Display script header
    echo "==================================================================="
    echo "SRA Download and Processing Script for L. pneumophila Analysis"
    echo "Version: 1.0"
    echo "Started at: $(date)"
    echo "==================================================================="
    
    # Store command line arguments
    SRA=$1
    REF_GEN_SITE=$2
    
    # Execute workflow steps
    setup_logging
    check_dependencies
    validate_inputs "$@"
    create_directories
    
    # Change to raw_data directory
    cd ./${SRA}/raw_data
    
    # Execute main processing steps
    echo "Step 1/4: Downloading SRA data..."
    download_sra
    
    echo "Step 2/4: Converting to FASTQ format..."
    convert_to_fastq
    
    echo "Step 3/4: Verifying file integrity..."
    verify_files
    
    echo "Step 4/4: Cleaning up temporary files..."
    cleanup
    
    # Return to original directory
    cd ../..
    
    # Display completion message
    echo "==================================================================="
    echo "Download process completed successfully at $(date)"
    echo "Output files can be found in ./${SRA}/raw_data/fastq/"
    echo "Log files are available in ./${SRA}/logs/"
    echo "==================================================================="
}

# Execute main function with command line arguments
main "$@"