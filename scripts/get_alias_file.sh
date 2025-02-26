#!/bin/bash

species="$1"
gca_id="$2"
output_dir="$3"
base_url_ncbi="https://hgdownload.soe.ucsc.edu/hubs"
success=false
gcf_used=false

# Function to download file using curl
download_file() {
    local url="$1"
    local output="$2"
    echo "$output"
    if curl -s -f -o "$output" "$url"; then
        # Check if the file contains HTML indicating a 404 error
        if grep -q "<!DOCTYPE HTML" "$output" && grep -q "<title>404 Not Found</title>" "$output"; then
            rm "$output"  # Remove the error page
            return 1
        else
            return 0
        fi
    else
        return 1
    fi
}

# Function to get corresponding GCF ID
get_gcf_id() {
    local gca_id="$1"
    local api_url="https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/${gca_id}"

    # Use jq to parse JSON response and extract paired_assembly_accession
    local gcf_id=$(curl -s "$api_url" | grep -o '"paired_assembly_accession":"GCF_[^"]*' | sed 's/"paired_assembly_accession":"//')
    echo "$gcf_id"
}

# Get GCF ID based on GCA ID
gcf_id=$(get_gcf_id "$gca_id")

# If GCF ID is empty, use GCA ID with 'A' replaced by 'F'
if [ -z "$gcf_id" ]; then
    gcf_id="${gca_id/A/F}"
fi

# Download GCA file first
part1=$(echo "$gca_id" | cut -c5-7)
part2=$(echo "$gca_id" | cut -c8-10)
part3=$(echo "$gca_id" | cut -c11-13)
url="${base_url_ncbi}/GCA/${part1}/${part2}/${part3}/${gca_id}/${gca_id}.chromAlias.txt"

if download_file "$url" "$output_dir/${gca_id}.chromAlias.txt"; then
    echo "Successfully downloaded GCA file for $species ($gca_id)"
    success=true
else
    # Attempt GCF download if GCA fails
    part1=$(echo "$gcf_id" | cut -c5-7)
    part2=$(echo "$gcf_id" | cut -c8-10)
    part3=$(echo "$gcf_id" | cut -c11-13)
    url="${base_url_ncbi}/GCF/${part1}/${part2}/${part3}/${gcf_id}/${gcf_id}.chromAlias.txt"
    if download_file "$url" "$output_dir/${gcf_id}.chromAlias.txt"; then
        echo "Successfully downloaded GCF file for $species ($gcf_id)"
        success=true
        gcf_used=true
    fi
fi

if ! $success; then
    echo "Failed to download for both GCA and GCF variants of $species"
    exit 1
fi

echo "Download process completed for $species."

