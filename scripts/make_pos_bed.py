#!/usr/bin/env python3

import polars as pl
import argparse
import os

def process_bed_file(input_file, outdir):
    # Read the BED file
    bed_data = pl.read_csv(input_file, separator='\t')

    # Create the 'hg38_position' range and drop 'start' and 'end'
    bed_data = bed_data.with_columns(
        pl.int_ranges('start', 'end').alias('hg38_position')
    ).drop('start', 'end')

    # Explode the position ranges
    bed_data = bed_data.explode('hg38_position')

    # Group by chromosome and position
    bed_data = bed_data.group_by('hg38_chr', 'hg38_position').all()

    # Join the 'region' column
    bed_data = bed_data.with_columns(
        pl.col('region').list.join(',').alias('region')
    )

    # Sort by chromosome and position
    bed_data = bed_data.sort(['hg38_chr', 'hg38_position'])

    # Adjust positions by adding 1
    bed_data = bed_data.with_columns(
        (pl.col('hg38_position') + 1).alias('hg38_position')
    )

    # Write to a Parquet file
    output_file = os.path.join(outdir, 'features_by_position.parquet')
    bed_data.write_parquet(output_file)

    print(f"Processed data saved to {output_file}")

def main():
    # Argument parser for input file
    parser = argparse.ArgumentParser(description="Process a BED file and generate a Parquet file by positions.")
    parser.add_argument("input_file", help="Path to the input BED file")
    parser.add_argument("outdir", help="Path to output directory")
    # Parse arguments
    args = parser.parse_args()

    # Process the BED file
    process_bed_file(args.input_file, args.outdir)

if __name__ == "__main__":
    main()

