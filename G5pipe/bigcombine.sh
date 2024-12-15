#!/bin/bash
# echo -e "SRA_FILE\tGROUP\tPOS\tREF\tALT\tQUAL\tGENE_NAME\tVARIANT_TYPE\tFEATURE_TYPE\tFEATURE_BIOTYPE" > $output_file
# Paths to inputs


# Paths to inputs
sragroup_file="SRAgroup.txt"
output_file="big_combine.vcf"

# Initialize the output file with the header
echo -e "SRA_RUN\tGROUP\tPOS\tREF\tALT\tQUAL\tGENE\tANNOTATION\tFEATURE_TYPE" > $output_file

# Iterate through each group listed in the SRAgroup.txt
while read -r group; do
    # Navigate to the group directory
    if [ -d "$group" ]; then
        echo "Processing group: $group"
        cd "$group" || { echo "Failed to enter directory $group"; continue; }

        # Find and process all annotated VCF files in the group directory
        find . -type f -name "annotated_*.vcf" | while read -r vcf_file; do
            # Extract the SRA_RUN from the filename
            SRA_RUN=$(basename "$vcf_file" | sed -E 's/annotated_(.+)\.vcf/\1/')

            # Process the VCF file using awk
            awk -v SRA_RUN="$SRA_RUN" -v group="$group" '
            BEGIN { OFS="\t" }
            /^##/ { next }  # Skip header lines
            /^#CHROM/ { next }  # Skip column header line
            {
                pos = $2
                ref = $4
                alt = $5
                qual = $6
                info = $8

                # Extract Gene_Name, Annotation, and Feature_Type from the ANN field in INFO
                match(info, /ANN=([^;]+)/, ann)
                split(ann[1], annotations, ",")

                gene_name = ""
                annotation = ""
                feature_type = ""

                for (i in annotations) {
                    split(annotations[i], fields, "|")
                    if (fields[4] != "") {
                        gene_name = fields[4]  # Gene_Name
                    }
                    if (fields[2] != "") {
                        annotation = fields[2]  # Annotation
                    }
                    if (fields[6] != "") {
                        feature_type = fields[6]  # Feature_Type
                    }

                    # Stop after finding the first valid set of values
                    if (gene_name != "" && annotation != "" && feature_type != "") {
                        break
                    }
                }

                # Print the selected fields with SRA_RUN and GROUP
                print SRA_RUN, group, pos, ref, alt, qual, gene_name, annotation, feature_type
            }
            ' "$vcf_file" >> "../$output_file"

        done

        # Return to the parent directory
        cd ..
    else
        echo "Group directory $group does not exist. Skipping."
    fi
done < "$sragroup_file"

# Sort the output file by POS column
sort -k3,3n $output_file -o $output_file

echo "Processing complete. Results saved to $output_file."
