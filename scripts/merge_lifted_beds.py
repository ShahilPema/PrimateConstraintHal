#!/usr/bin/env python3

import polars as pl
import sys
import os
from pathlib import Path
from joblib import Parallel, delayed

def process_and_save_batch(bed_files, species, batch_index):
    # Lazily read and concatenate all BED files in the batch
    lifted_data = pl.scan_csv(bed_files, separator='\t', has_header=False, raise_if_empty=False)

    # Process the data: rename columns and extract information
    lifted_data = lifted_data.with_columns(
        (pl.col('column_1')).alias('hg38_chr'),
        pl.col('column_3').alias('hg38_position'),
        pl.col('column_4').str.split(':').list.get(0).alias(f'{species}_contig'),
        pl.col('column_4').str.split(':').list.get(1).alias(f'{species}_position')
    ).select(['hg38_chr', 'hg38_position', f'{species}_contig', f'{species}_position'])

    # Collect and write to a Parquet file
    batch_file = f'batch_{batch_index}_{species}_alignment.parquet'
    lifted_data.sink_parquet(batch_file)
    print(f"Processed batch {batch_index} written to {batch_file}")

def merge_parquet_files(species):
    # Gather all Parquet files and merge them
    parquet_files = sorted(Path('.').glob(f'batch_*_{species}_alignment.parquet'))
    if not parquet_files:
        print("No Parquet files found to merge.")
        return

    merged_df = pl.scan_parquet(parquet_files)
    merged_output_file = f'{species}_alignment.parquet'
    merged_df.sink_parquet(merged_output_file)
    print(f"All batches merged into {merged_output_file}")

    # Cleanup individual batch files
    for parquet_file in parquet_files:
        #os.remove(parquet_file)
        print(f"Deleted batch file {parquet_file}")

def main():
    if len(sys.argv) < 3:
        print("Usage: ./merge_lifted_beds.py <species> <bed_file1> <bed_file2> ...")
        sys.exit(1)
    
    species = sys.argv[1]
    bed_files = sys.argv[2:]
    
    # Process files in batches of 200
    batch_size = 200
    Parallel(n_jobs=-2)(delayed(process_and_save_batch)(bed_files[i:i+batch_size], species, i // batch_size) for i in range(0, len(bed_files), batch_size))

    # Merge all the batch Parquet files
    merge_parquet_files(species)

if __name__ == "__main__":
    main()

