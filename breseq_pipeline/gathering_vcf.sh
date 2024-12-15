#!/bin/bash
# This script gathers all the output vcf file to a directory
# Author: Jason
# Date: 25/11/2024

# Define Args
RESULTS_DIR="/hdd1/home/f24_wyli/project/pipeline/analysis_results/results" 
OUTPUT_DIR="/hdd1/home/f24_wyli/project/pipeline/vcf_gathering"   
# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Iterate over all child dirs of SRR parent dir
for SRR_DIR in "${RESULTS_DIR}"/SRR*; do
    if [ -d "${SRR_DIR}" ]; then
        # Get SRR Number
        SRR=$(basename "${SRR_DIR}")
        
        # Look for all vcf file
        VCF_FILE=$(find "${SRR_DIR}/output/output.vcf" -maxdepth 1 -type f -name "*.vcf")
        
        if [ -f "${VCF_FILE}" ]; then
            # copy and rename
            cp "${VCF_FILE}" "${OUTPUT_DIR}/${SRR}.vcf"
            echo "Copy ${VCF_FILE} to ${OUTPUT_DIR}/${SRR}.vcf"
        else
            echo "Not Found VCF IN ${SRR_DIR} "
        fi
    fi
done

