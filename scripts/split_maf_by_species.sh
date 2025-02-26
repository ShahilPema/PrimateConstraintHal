#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_maf> <species> <output_file>"
    exit 1
fi

# Assign the arguments to variables for better readability
input_maf="$1"
species="$2"
output_file="$3"

# Create a temporary species file with the species name
temp_species_file=$(mktemp)
echo "hg38" > "$temp_species_file"
echo "$species" >> "$temp_species_file"

# Run the mafSpeciesSubset command
mafSpeciesSubset "$input_maf" "$temp_species_file" "$output_file"

# Clean up the temporary file
rm "$temp_species_file"
