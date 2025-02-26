import polars as pl
import sys
import os

# Inputs from the command line
species = sys.argv[1]
lifttable_path = sys.argv[2]
sample_parquet_path = sys.argv[3]
output_dir = sys.argv[4]

# Extract sampleid from the filename (assuming it's the part before .parquet)
sampleid = sample_parquet_path.split('/')[-1].replace('_filtered.parquet', '')

# Read the sample Parquet file and add a count column specific to the sampleid
sample = pl.scan_parquet(sample_parquet_path).with_columns(
        pl.col('ALT').str.to_uppercase().alias(f"{species}_ALT"),
        ).drop('ALT','REF')
lifttable = pl.scan_parquet(lifttable_path).with_columns(
        pl.col('ALT').str.to_uppercase().alias("hg38_ALT"),
        pl.col('REF').str.to_uppercase().alias("hg38_REF"),
        pl.col(f"{species}_REF").str.to_uppercase().alias(f"{species}_REF")
        ).drop('ALT','REF')

counts = lifttable.join(sample, how='left', on=[f"{species}_contig", f"{species}_position"], coalesce=True)
counts = counts.with_columns(
    pl.when(
        (pl.col(f"{species}_ALT").is_not_null()) &
        (pl.col(f"{species}_ALT") == pl.col("hg38_ALT")) & 
        (pl.col(f"{species}_REF") == pl.col("hg38_REF"))
         )
    .then(pl.lit(1))
    .when(
        (pl.col(f"{species}_ALT").is_not_null()) &
        (pl.col(f"{species}_ALT") != pl.col("hg38_ALT")) &
        (pl.col(f"{species}_REF") == pl.col("hg38_REF"))
        )
    .then(pl.lit(0))
    .when(
        (pl.col(f"{species}_ALT").is_null()) &
        (pl.col(f"{species}_contig").is_not_null())
        )
    .then(pl.lit(0))
    .when(
        (pl.col(f"{species}_ALT").is_not_null()) &
        (pl.col(f"{species}_REF") != pl.col("hg38_REF"))
        )
    .then(pl.lit(2))
    .otherwise(pl.lit(None))
    .alias(f"{species}_{sampleid}")
    )

counts = counts.drop(f"{species}_contig", f"{species}_position", f"{species}_REF", f"{species}_ALT", "hg38_ALT", "hg38_REF")

# Write the final output file
output_file = f"{output_dir}/{sampleid}_counts.parquet"
counts.sink_parquet(output_file)
