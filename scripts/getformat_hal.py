#!/usr/bin/env python3

import subprocess
import polars as pl
import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: script.py <species> <alias_file> <hal_file>")
        sys.exit(1)

    # Parse command-line arguments
    species = sys.argv[1]
    alias_file = sys.argv[2]
    hal_file = sys.argv[3]

    # Load alias dictionary (modify as necessary to load your alias data)
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

    # Run HAL command to get contigs
    hal_command = f"halStats --sequences {species} {hal_file}"
    result = subprocess.run(hal_command, shell=True, capture_output=True, text=True, check=True)
    hal_contigs = result.stdout.split(',')

    # Process HAL contigs and determine the most frequent assembly
    hal_data = pl.DataFrame(hal_contigs, schema=pl.Schema({'column_1': pl.Utf8}))
    hal_assembly = hal_data.select(pl.col('column_1')).to_series().replace_strict(alias_dict, default=None).value_counts().sort('count', descending=True).select('column_1').to_series().item(0)

    print(hal_assembly)

def load_alias_dict(alias_file):
    # Placeholder for loading alias dictionary
    # Implement logic to parse the alias file into a dictionary
    return {}

if __name__ == "__main__":
    main()
