#!/bin/bash
species=$1
hal=$2

# Extract FASTA
hal2fasta ${hal} ${species} --outFaPath "${species}.fa"

# Create FAI index
samtools faidx "${species}.fa"

