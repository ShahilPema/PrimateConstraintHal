import polars as pl
import sys
import os

def merge_parquets(species, parquet_files, output_dir):
    # Ensure there are parquet files to process
    if not parquet_files:
        print(f"No parquet files provided for species {species}")
        return
    
    # Read and stack all parquet files
    dfs = [pl.read_parquet(file) for file in parquet_files]
    df = pl.concat(dfs, how="vertical")  # Vertically stack the DataFrames
    df = df.unique()

    # Generate species-specific output filename and path
    output_filename = f"{species}_alignment.parquet"
    output_path = os.path.join(output_dir, output_filename)
    
    # Write the final DataFrame to a parquet file
    df.write_parquet(output_path)
    
    print(f"Processed {len(parquet_files)} files for species: {species}")
    print(f"Final dataframe shape: {df.shape}")

if __name__ == "__main__":
    species = sys.argv[1]
    parquet_files = sys.argv[2].split(',')
    output_dir = sys.argv[3]  # Output directory passed as third argument
    merge_parquets(species, parquet_files, output_dir)
