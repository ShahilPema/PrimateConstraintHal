import polars as pl
import argparse
import glob
import os

def merge_bed_files(bed_files, pos_bed, sample_id,  output_prefix):
    # Read all BED files in parallel, assuming no headers in the files
    dfs = [pl.read_csv(file, has_header=False, separator="\t") for file in bed_files if os.path.getsize(file) > 0]

    # Concatenate all DataFrames into one
    merged_df = pl.concat(dfs)
    
    # Save the first few rows to a CSV file
    output_csv = f"{output_prefix}_head.csv"
    merged_df.head(10).write_csv(output_csv, include_header=False, separator="\t")
    print(f"First 10 rows written to {output_csv}")

    positions = pl.read_parquet(pos_bed).select('hg38_chr', 'hg38_position')
    
    merged_df = merged_df.select('column_1', 'column_3').unique().with_columns(
        pl.lit(1).alias('variant_binary')
    )
    merged_df.columns = ['hg38_chr', 'hg38_position', sample_id]
    
    positions = positions.join(merged_df, on=['hg38_chr', 'hg38_position'], coalesce=True, how='left')
 
    # Write the merged DataFrame to a Parquet file
    output_parquet = f"{output_prefix}.parquet"
    merged_df.write_parquet(output_parquet)
    print(f"Merged data written to {output_parquet}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Merge BED files using Polars.")
    parser.add_argument("bed_files", nargs='+', help="List of BED files to merge")
    parser.add_argument("-p", "--posbed", required=True, help="Position parquet to merge to")
    parser.add_argument("-i", "--sampleid", required=True, help="ID for the column name")
    parser.add_argument("-o", "--output", required=True, help="Output prefix for the merged files")
    
    args = parser.parse_args()
    
    # Merge the BED files and write outputs
    merge_bed_files(args.bed_files, args.posbed,  args.sampleid, args.output)

if __name__ == "__main__":
    main()

