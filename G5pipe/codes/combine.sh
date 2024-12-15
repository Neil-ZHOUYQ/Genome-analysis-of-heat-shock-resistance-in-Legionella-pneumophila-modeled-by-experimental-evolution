#!/bin/bash

# Define the root directory containing subdirectories (SRR_1, SRR_2, etc.)
root_dir=$1
output_file=$2

group=$(basename "$(pwd)")



# Append the new column header for the output
echo -e "IDENTIFIER\tGROUP\tPOS\tREF\tALT\tQUAL\tGENE" > $output_file

# Process each VCF file for data rows
find "$root_dir" -type f -name "annotated_*.vcf" | while read -r vcf_file; do
    # Extract the IDENTIFIER from the filename
    identifier=$(basename "$vcf_file" | sed -E 's/annotated_(.+)\.vcf/\1/')

    awk -v identifier="$identifier" -v group="$group" '
    BEGIN { OFS="\t" }
    /^##/ { next }  # Skip header lines
    /^#CHROM/ { next }  # Skip column header line
    {
        pos = $2
        ref = $4
        alt = $5
        qual = $6
        info = $8

        # Extract Gene_Name from the ANN field in INFO
        match(info, /ANN=([^;]+)/, ann)
        split(ann[1], annotations, ",")
        gene_name = ""
        for (i in annotations) {
            split(annotations[i], fields, "|")
            if (fields[4] != "") {
                gene_name = fields[4]
                break
            }
        }

        # Print the selected fields with IDENTIFIER and GROUP
        print identifier, group, pos, ref, alt, qual, gene_name
    }
    ' "$vcf_file" >> $output_file

done

# Sort the output file by POS column
sort -k3,3n $output_file -o $output_file

echo "Processing complete. Results saved to $output_file."
