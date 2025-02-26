import polars as pl
import argparse
import os
from datetime import datetime

def merge_species_counts(output_file, count_files):
    batch_size = 50
    temp_files = []

    # Process all files except the first one in batches
    remaining_files = count_files[1:]
    for i in range(0, len(remaining_files), batch_size):
        batch = remaining_files[i:i+batch_size]
        
        temp_file_name = f"tmp{len(temp_files)+1}.parquet"
        #dfs = [pl.read_parquet(file).select(pl.nth(2)) for file in batch]
        
        #merged_batch = pl.concat(dfs, how="horizontal", rechunk=True)
        #merged_batch.write_parquet(temp_file_name)
        temp_files.append(temp_file_name)

    # Final merge including the first file
    first = pl.scan_parquet(count_files[0])
    final_dfs = [first] + [pl.scan_parquet(file) for file in temp_files]
    final_dfs = final_dfs[0:2]
    final_merge = pl.concat(final_dfs, how="horizontal")
    final_merge.collect().write_parquet(output_file)

    # Clean up temporary files
    #for temp_file in temp_files:
    #    os.remove(temp_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge species count files into a single Parquet file.")
    parser.add_argument('count_files', help='Comma-separated list of count files to merge')
    
    args = parser.parse_args()
    
    # Get current timestamp
    timestamp = datetime.now().strftime("%Y%m%d")
    output_file = f"merged_counts_{timestamp}.parquet"

    # Split count files into a list
    count_files_list = args.count_files.split(',')
    
    merge_species_counts(output_file, count_files_list)

