#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <bed_file> <species> <hal_file> -o <output_directory>"
    exit 1
fi

# Assign input arguments to variables
bed_file="$1"
species="$2"
hal_file="$3"

# Parse the -o flag and output directory
if [ "$4" == "-o" ]; then
    output_dir="$5"
else
    echo "Error: Missing or incorrect -o flag"
    exit 1
fi

# Check if the input files exist
if [ ! -f "$bed_file" ]; then
    echo "Error: BED file $bed_file not found."
    exit 1
fi

if [ ! -f "$hal_file" ]; then
    echo "Error: HAL file $hal_file not found."
    exit 1
fi

# Check if the output directory exists, if not, create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Extract basename of the BED file (without extension)
bed_basename=$(basename "$bed_file" .bed)

# Step 1: Add a 4th column to the BED file with concatenated coordinates
awk 'BEGIN {OFS="\t"} {print $0, $1":"$3}' "$bed_file" > "${output_dir}/${bed_basename}_annotated.bed"

# Step 2: Lift BED from human to species
if ! halLiftover "$hal_file" Homo_sapiens "${output_dir}/${bed_basename}_annotated.bed" "$species" "${output_dir}/${bed_basename}_species_lifted.bed" --bedType 3; then
    echo "Error: halLiftover failed during human to species liftover."
    exit 1
fi

# Step 3: Add a 5th column to the lifted BED file with species coordinates
awk 'BEGIN {OFS="\t"} {print $0, $1":"$3}' "${output_dir}/${bed_basename}_species_lifted.bed" > "${output_dir}/${bed_basename}_species_lifted_annotated.bed"

# Step 4: Lift the annotated BED back to human
if ! halLiftover "$hal_file" "$species" "${output_dir}/${bed_basename}_species_lifted_annotated.bed" Homo_sapiens "${output_dir}/${bed_basename}_human_lifted.bed" --bedType 3; then
    echo "Error: halLiftover failed during species to human liftover."
    exit 1
fi

# Step 5: Add a 5th column to the human lifted BED file with human coordinates
awk 'BEGIN {OFS="\t"} {print $4, $5, $1":"$3}' "${output_dir}/${bed_basename}_human_lifted.bed" > "${output_dir}/${bed_basename}_final.bed"

# Cleanup intermediate files
rm -f "${output_dir}/${bed_basename}_annotated.bed" \
     "${output_dir}/${bed_basename}_species_lifted.bed" \
     "${output_dir}/${bed_basename}_species_lifted_annotated.bed" \
     "${output_dir}/${bed_basename}_human_lifted.bed"