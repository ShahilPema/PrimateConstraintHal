#!/bin/bash

# Run bs list dataset and extract the ID for primate_genomic_variation
dataset_id=$(bs list dataset --project-id 393739348 | awk '/primate_genomic_variation/ {print $4}')

# Check if dataset_id was found
if [ -z "$dataset_id" ]; then
    echo "Error: Could not find dataset ID for primate_genomic_variation"
    exit 1
fi

# Run bs dataset contents with the extracted ID and save the output to a file
bs dataset contents -i ds.13ff18d1208042048fd3066fd6cf673d | grep ".vcf.gz" | grep -v ".tbi" | awk '{print $2 "\t" $4}' | sed 's/variant_calls_SNV_v1\///g' > vcf_ids.txt

count=$(wc -l < vcf_ids.txt)
echo "Extracted $count file IDs"
