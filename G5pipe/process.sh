#!/bin/bash

# Paths
SRAgroup_file="./SRAgroup.txt"
SRArun_file="./sorted_SRArun.txt"
source_path="./reference"  # Path where the directories to copy are located
current_path=$(pwd)

# Read group names
group_names=$(cat "$SRAgroup_file")

# Create directories for each group and process
while read -r group; do
    # Create the group directory
    mkdir -p "$current_path/$group"

    # Copy the required directories into the group directory
    cp  -r ./codes "$current_path/$group/"
    cp  -r ./reference "$current_path/$group/"

    # Find all SRA identifiers for this group
    grep -P "\t$group$" "$SRArun_file" | cut -f1 | while read -r sra_id; do
        # Run the command script with the SRA identifier
        (cd "$current_path/$group"
        ./codes/commands.sh "$sra_id"
        )
    done
    (cd "$current_path/$group"
    ./codes/combine.sh . ./combined.vcf
    )
    
done <<< "$group_names"

echo "Whole processing complete."
