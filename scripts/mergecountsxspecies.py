import polars as pl
import os
import argparse
from pathlib import Path

def process_batch(species, batch_idx, parquet_files, output_dir):
    """
    Processes a batch of count files and writes the concatenated result as a Parquet file.

    Args:
        species (str): Species name.
        batch_idx (int): Batch index.
        parquet_files (list): List of Parquet file paths for this batch.
        output_dir (str): Directory to save the intermediate output files.
    """
    if batch_idx == 0:
        first_df = pl.scan_parquet(parquet_files[0])
        dfs = [pl.scan_parquet(file).select(pl.nth(-1)) for file in parquet_files[1:]]
        dfs = [first_df] + dfs
    else:
        dfs = [pl.scan_parquet(file).select(pl.nth(-1)) for file in parquet_files]

    concatenated_df = pl.concat(dfs, how="horizontal")
    batch_output = os.path.join(output_dir, f"{species}_batch_{batch_idx}.parquet")
    concatenated_df.collect().write_parquet(batch_output)


def merge_batches(species, num_batches, output_dir):
    """
    Merges all batch files into a single Parquet file.

    Args:
        species (str): Species name.
        num_batches (int): Number of batches to merge.
        output_dir (str): Directory containing batch files.
    """
    batch_files = [os.path.join(output_dir, f"{species}_batch_{i}.parquet") for i in range(num_batches)]
    dfs = [pl.scan_parquet(file) for file in batch_files]
    merged_df = pl.concat(dfs, how="horizontal")
    output_filename = os.path.join(output_dir, f"{species}_counts_temp.parquet")
    merged_df.collect().write_parquet(output_filename)
    print(f"Merged {num_batches} batches into: {output_filename}")


def reposition_counts(species, output_dir, bed_file_path):
    """
    Repositions counts based on BED file data, partitioned by a specified column.

    Args:
        species (str): Species name.
        output_dir (str): Directory containing the merged counts file.
        bed_file_path (str): Path to the BED position file.
        partition_col (str): Column to partition the output file by (e.g., chromosome).
    """
    counts_file = os.path.join(output_dir, f"{species}_counts")

    # Read data lazily
    counts_data = pl.scan_parquet(f"{counts_file}_temp.parquet")
    bed_data = pl.scan_parquet(bed_file_path)

    # Perform the join
    joined_data = bed_data.join(
        counts_data,
        on=["hg38_chr", "hg38_position", "ALT"],
        how="left",
    )

    # Write output partitioned by the specified column
    joined_data.sink_parquet(f"{counts_file}.parquet")

def reposition_counts(species, output_dir, bed_file_path):
    """
    Repositions counts based on BED file data using lazy evaluation and batch processing.

    Args:
        species (str): Species name.
        output_dir (str): Directory containing the merged counts file.
        bed_file_path (str): Path to the BED position file.
    """
    counts_file = os.path.join(output_dir, f"{species}_counts")
    batch_dir = os.path.join(output_dir, "batches")
    Path(batch_dir).mkdir(exist_ok=True)

    # Scan bed data to get total number of rows
    bed_data = pl.scan_parquet(bed_file_path)
    total_rows = bed_data.select(pl.len()).collect().item()

    # Process in batches of 1 million rows
    batch_size = 100_000_000
    for batch_num, offset in enumerate(range(0, total_rows, batch_size)):
        print(f"Processing batch {batch_num + 1}")

        # Get current batch of bed data lazily
        batch_bed = bed_data.slice(offset, batch_size)

        # Get unique chromosome and position combinations for this batch
        chroms = batch_bed.select(pl.col("hg38_chr").unique()).collect().to_series().to_list()
        max_pos = batch_bed.select(pl.col("hg38_position").max()).collect().item()
        min_pos = batch_bed.select(pl.col("hg38_position").min()).collect().item()
        
        # Scan and filter counts data for just this batch
        batch_counts = (pl.scan_parquet(f"{counts_file}_temp.parquet")
            .filter(
                (pl.col("hg38_chr").is_in(chroms)),
                (pl.col("hg38_position").is_between(min_pos, max_pos))
                )
        )

        # Perform lazy join for this batch
        batch_joined = batch_bed.join(
            batch_counts,
            on=["hg38_chr", "hg38_position", "ALT"],
            how="left"
        )

        # Write batch result
        batch_file = os.path.join(batch_dir, f"batch_{batch_num:04d}.parquet")
        batch_joined.sink_parquet(batch_file)

    # Combine all batches lazily
    print("Combining batches...")
    combined_data = pl.scan_parquet(f"{batch_dir}/*.parquet")
    combined_data.sink_parquet(f"{counts_file}.parquet")

    # Clean up batch files
    if os.path.exists(f"{counts_file}.parquet"):
        import shutil
        shutil.rmtree(batch_dir)
        print("Batch files cleaned up")
    else:
        print("Warning: Final file not created, keeping batch files for inspection")

def process_species(species, count_files, bed_file, output_dir, batch_size=10):
    """
    Processes counts for a species.

    Args:
        species (str): Species name.
        count_files (str): Comma-separated list of count file paths.
        bed_file (str): Path to the BED position file.
        output_dir (str): Directory to save the output Parquet files.
        batch_size (int): Number of files to process per batch.
    """
    os.makedirs(output_dir, exist_ok=True)
    parquet_files = [file.strip() for file in count_files.split(",")]

    if not parquet_files:
        print(f"No count files provided for species: {species}")
        return

    num_batches = (len(parquet_files) + batch_size - 1) // batch_size
    for batch_idx in range(num_batches):
        batch_files = parquet_files[batch_idx * batch_size:(batch_idx + 1) * batch_size]
        #process_batch(species, batch_idx, batch_files, output_dir)
        print(f"Processed batch {batch_idx}")

    #merge_batches(species, num_batches, output_dir)
    reposition_counts(species, output_dir, bed_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge count files by species.")
    parser.add_argument("species", help="The species name.")
    parser.add_argument("count_files", help="Comma-separated list of count files.")
    parser.add_argument("bed_file", help="Path to the BED position file.")
    parser.add_argument("output_dir", help="Directory to save the output Parquet files.")

    args = parser.parse_args()
    process_species(args.species, args.count_files, args.bed_file, args.output_dir)

