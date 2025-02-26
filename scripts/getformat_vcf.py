#!/usr/bin/env python3

import subprocess
import os
import polars as pl
import sys

# Ensure the correct number of arguments
if len(sys.argv) != 3:
    print("Usage: process_vcf.py <alias_file> <vcf_file>")
    sys.exit(1)

vcf_file = sys.argv[1]
alias_file = sys.argv[2]

# Step 1: Load the alias file and create a dictionary
try:
    alias_data = pl.read_csv(alias_file, separator='\t')
    alias_data.columns = [item.replace(" ", "").replace("#", "") for item in alias_data.columns]
    alias_dict = {}
    for column in alias_data.columns:
        contig_format = alias_data.select(column)
        contig_format = contig_format.with_columns(
            pl.lit(column).alias('format')
        )
        alias_dict.update(dict(contig_format.iter_rows()))
except FileNotFoundError:
    print(f"Alias file not found: {alias_file}")
    sys.exit(1)
except Exception as e:
    print(f"Unexpected error occurred while creating alias dictionary: {e}")
    sys.exit(1)

# Step 2: Determine assemblies for HAL and VCF
try:
    vcf_data = pl.read_csv(vcf_file, separator='\t', comment_prefix='#', has_header=False)
    vcf_assembly = (
        vcf_data
        .select(pl.col('column_1'))
        .unique()
        .to_series()
        .replace_strict(alias_dict, default=None)
        .value_counts()
        .sort('count', descending=True)
        .select('column_1')
        .to_series()
        .item(0)
    )
    print(vcf_assembly if vcf_assembly is not None else 'Other')
except FileNotFoundError:
    print(f"VCF file not found: {vcf_file}")
    sys.exit(1)
except Exception as e:
    print(f"Error determining assemblies for HAL and VCF: {e}")
    sys.exit(1)

