#!/usr/bin/env python3

import sys
import polars as pl

def main(alignment_file, positions_file):
    # Load positions
    positions = pl.scan_parquet(positions_file)
    
    # Extract unique hg38 chromosomes from alignment
    chrs = pl.scan_parquet(alignment_file).select('hg38_chr').unique().collect().to_series().to_list()

    # Collect data for each chromosome
    all_data = []
    for chrm in chrs:
        df = pl.scan_parquet(alignment_file).filter(pl.col('hg38_chr') == chrm).unique()
        all_data.append(df)

    # Concatenate the data for all chromosomes
    alignment = pl.concat(all_data, how='vertical')

    # Perform the join on hg38_chr and hg38_position
    merged = positions.join(alignment, on=['hg38_chr', 'hg38_position'], how='left', coalesce=True)

    # Save the result
    output_file = f"{alignment_file.split('.')[0]}_v2.parquet"
    merged.collect().write_parquet(output_file)
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: add_all_pos.py <alignment_file> <positions_file>")
        sys.exit(1)
    alignment_file = sys.argv[1]
    positions_file = sys.argv[2]

    main(alignment_file, positions_file)
