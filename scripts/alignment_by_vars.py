import polars as pl
import sys
from pyfaidx import Fasta

# Inputs from the command line
species = sys.argv[1]
lifttable_path = sys.argv[2]
vars_parquet_path = sys.argv[3]
fasta_path = sys.argv[4]
output_dir = sys.argv[5]

fasta = Fasta(fasta_path)

lifttable = pl.read_parquet(lifttable_path)

lifttable = lifttable.with_columns(
        pl.col(f"{species}_position").cast(pl.Int64).alias(f"{species}_position"),
        pl.col('hg38_position').cast(pl.Int64).alias('hg38_position')
)

lifttable = lifttable.filter(lifttable.select("hg38_chr", "hg38_position").is_duplicated().not_())

# Add the reference allele to the dataframe
lifttable = lifttable.with_columns(
    pl.struct([
        pl.col(f"{species}_contig").alias("contig"),
        (pl.col(f"{species}_position").cast(pl.Int64) - 1).alias("position")
    ]).map_elements(lambda x: fasta[x["contig"]][x["position"]].seq, return_dtype=pl.Utf8)
    .alias(f"{species}_REF")
)

vars_table = pl.read_parquet(vars_parquet_path)
lifttable = vars_table.join(lifttable, how='left', on=["hg38_chr", "hg38_position"], coalesce=True)

lifttable = lifttable.with_columns(
        pl.when(pl.col(f"{species}_REF") != pl.col("REF"))
        .then(pl.lit(None))
        .otherwise(pl.col("ALT"))
        .alias("ALT")
        )

lifttable.write_parquet(output_dir)
