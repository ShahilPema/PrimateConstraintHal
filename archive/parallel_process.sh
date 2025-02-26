#!/bin/bash

# List of files to search for
files=(
    PD_0402.SNV.vcf.gz PD_0078.SNV.vcf.gz PD_0358.SNV.vcf.gz
    PD_0371.SNV.vcf.gz PD_0361.SNV.vcf.gz PD_0311.SNV.vcf.gz
    PD_0012.SNV.vcf.gz PD_0373.SNV.vcf.gz PD_0374.SNV.vcf.gz
    PD_0011.SNV.vcf.gz PD_0375.SNV.vcf.gz PD_0353.SNV.vcf.gz
    PD_0370.SNV.vcf.gz PD_0369.SNV.vcf.gz PD_0354.SNV.vcf.gz
    PD_0372.SNV.vcf.gz PD_0360.SNV.vcf.gz PD_0013.SNV.vcf.gz
    PD_0079.SNV.vcf.gz PD_0077.SNV.vcf.gz PD_0007.SNV.vcf.gz
    PD_0141.SNV.vcf.gz PD_0363.SNV.vcf.gz PD_0362.SNV.vcf.gz
    PD_0366.SNV.vcf.gz PD_0359.SNV.vcf.gz PD_0357.SNV.vcf.gz
    PD_0367.SNV.vcf.gz PD_0355.SNV.vcf.gz PD_0356.SNV.vcf.gz
    PD_0310.SNV.vcf.gz PD_0368.SNV.vcf.gz
)

# Base directory to search
base_dir="./work/*/*/"

# Create a function to process each file
process_file() {
    local file="$1"
    found_file=$(find . -type f -name "${base_dir}$file" ! -xtype l)

    if [[ -n $found_file ]]; then
        # Count non-header lines
        count=$(zcat "$found_file" | grep -v '^#' | wc -l)
        echo "Non-header line count for $file: $count"
    else
        echo "$file not found!"
    fi
}

export -f process_file

# Send the file list to xargs for parallel processing with 4 processes
printf "%s\n" "${files[@]}" | xargs -P 4 -I {} bash -c 'process_file "$@"' _ {}

