#!/bin/bash

# Input parameters
CHR=$1
BED_FILE=$2
INPUT_MAF=$3
REF_SPECIES_FILE=$4


# Set file names for this chromosome
BED_NAME=$(basename "$BED_FILE" .bed) 
TMP_BED_FILE="chr${CHR}_tmp.bed"
TMP_BED_POS="chr${CHR}_tmp_bed_pos.txt"
TMP_MAF_POS="chr${CHR}_tmp_maf_pos.txt"
COMP_FILE="chr${CHR}_comp.txt"
REGIONS_MAF="chr${CHR}_regions.maf"
FINAL_MAF="${BED_NAME}_extracted.maf"

# Create temporary bed file with end positions incremented by 1
awk 'BEGIN {OFS="\t"} {print $1, $2+1, $3+1}' "${BED_FILE}" > "${TMP_BED_FILE}"

# Run mafsInRegion with the temporary bed file
mafsInRegion "${TMP_BED_FILE}" "${REGIONS_MAF}" "${INPUT_MAF}"

# Run mafSpeciesSubset
mafSpeciesSubset "${REGIONS_MAF}" "${REF_SPECIES_FILE}" "${FINAL_MAF}"

# Extract Positions from Bed File
awk '{for (i=$2; i<$3; i++) print $1":"i}' "${TMP_BED_FILE}" | sort -u > "${TMP_BED_POS}"

# Extract positions from MAF file
awk '/^s hg38.chr${CHR}/ {
    start=$3; 
    end=start+$4-1; 
    for (i=start; i<=end; i++) print "chr${CHR}:"i
}' "${FINAL_MAF}" | sort -u > "${TMP_MAF_POS}"

# Compare bed positions with maf positions
comm -23 "${TMP_BED_POS}" "${TMP_MAF_POS}" > "${COMP_FILE}"

line_count=$(wc -l < "${COMP_FILE}")

if [ "$line_count" -eq 0 ]; then
    echo "All positions in bed file were found in maf for chromosome ${CHR}"
    rm "${COMP_FILE}"
else
    echo "ERROR: $line_count positions in bed file were not in maf for chromosome ${CHR}"
fi

# Remove the temporary files
#rm "${TMP_BED_FILE}" "${TMP_BED_POS}" "${TMP_MAF_POS}" 
