#!/usr/bin/env python3
import sys
import os
#os.environ["POLARS_MAX_THREADS"] = "2"
import polars as pl
import gc
from pathlib import Path

def info2dict(info_string):
    data_dict = dict(item.split('=') for item in info_string.split(';'))
    default_keys = ['AC', 'AF', 'AN', 'BaseQRankSum', 'DP', 'ExcessHet', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR']
    data_list = []
    for key in default_keys:
        data_list.append(data_dict.get(key, None))
    return data_list

def filter_parquet(species, input_file):
    input_path = Path(input_file)

    output_file = input_path.with_name(f"{input_path.stem}_filtered.parquet")

    (pl.scan_parquet(input_path)
    .with_columns(
        pl.col('SAMPLE').str.splitn(":",3).struct.rename_fields(['GT', 'AD', 'OTHER']).alias('SAMPLE_STRUCT'),
        pl.col('INFO').map_elements(info2dict, return_dtype=pl.List(pl.String))
                       .list.to_struct(fields=['AC', 'AF', 'AN', 'BaseQRankSum', 'DP', 'ExcessHet', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR'])
                       .alias("INFO_STRUCT")
    )
    .unnest(['INFO_STRUCT', 'SAMPLE_STRUCT'])
    .drop('OTHER')
    .filter(
        (
            (pl.col('QD').cast(pl.Float64) >= 2) &
            (pl.col('FS').cast(pl.Float64) <= 60) &
            (pl.col('MQ').cast(pl.Float64) >= 40) &
            (pl.col('SOR').cast(pl.Float64) <= 3)
        ) & (
            (
                ((pl.col('GT') == "0/1") | (pl.col('GT') == "0|1")) &
                (pl.col('ReadPosRankSum').cast(pl.Float64) >= -8) &
                (pl.col('MQRankSum').cast(pl.Float64) >= -12.5) &
                (pl.col('AD').str.split(",").list.get(1).cast(pl.Int64) > 3)
            ) | (
                (pl.col('GT') == "1/1") | (pl.col('GT') == "1|1")
            )
        )
    )
    .select([pl.col(f"{species}_contig"), pl.col(f"{species}_position"), pl.col('REF'), pl.col('ALT')])
    .sink_parquet(output_file)
    )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: filter_script.py <species> <input_parquet_file>")
        sys.exit(1)

    species = sys.argv[1]
    input_file = sys.argv[2]
    filter_parquet(species, input_file)

    gc.collect()

