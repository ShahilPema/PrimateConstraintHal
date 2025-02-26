#!/usr/bin/env python3

import polars as pl
import sys
import os
from joblib import Parallel, delayed

def process_chromosome(chr_df, output_prefix, chr_name, max_rows=100000):
    num_partitions = (len(chr_df) + max_rows - 1) // max_rows
    for i in range(num_partitions):
        start_idx = i * max_rows
        length = min(max_rows, len(chr_df) - start_idx)
        # Slice the dataframe
        partition = chr_df.slice(start_idx, length)

        # Write partition to a BED file
        output_file = f"{output_prefix}_{chr_name}_part_{i + 1}.bed"
        partition.write_csv(output_file, include_header=False, separator="\t")
        print(f"Partition {i + 1} for chromosome {chr_name} written to {output_file}")

def split_and_write_bed_by_chr(df, output_prefix, max_rows=100000):
    # Partition the dataframe by chromosome using Polars' partition_by
    chr_partitions = df.partition_by("hg38_chr")

    # Use Parallel to process each chromosome in parallel with n_jobs set to the number of chromosomes
    Parallel(n_jobs=len(chr_partitions))(
        delayed(process_chromosome)(
            chr_df,
            output_prefix,
            chr_df['hg38_chr'][0],  # Get the chromosome name
            max_rows
        ) for chr_df in chr_partitions
    )

def main():
    if len(sys.argv) != 3:
        print("Usage: ./split_parquet_to_bed.py <input_parquet_file> <output_prefix>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]

    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} not found.")
        sys.exit(1)

    # Read the Parquet file
    bed_data = pl.read_parquet(input_file)
    bed_data = bed_data.with_columns(
        (pl.col('hg38_position') - 1).alias('start'),
        pl.col('hg38_position').alias('end')
    )
    bed_data = bed_data.select('hg38_chr', 'start', 'end')

    # Split and write the dataframe into BED files by chromosome
    split_and_write_bed_by_chr(bed_data, output_prefix)

if __name__ == "__main__":
    main()

