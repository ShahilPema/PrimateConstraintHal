import hail as hl
import argparse
import psutil
import math

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process a VCF file with Hail.")
parser.add_argument("--species", required=True, help="Species")
parser.add_argument("--cpus", required=True, help="Total cpus on machine")
parser.add_argument("--mt_files", required=True, help="List of mt files to merge")
parser.add_argument("--output", required=True, help="Output directory")
args = parser.parse_args()

hl_cpus = 64
forks = math.floor(int(args.cpus)/hl_cpus)

memory = int(psutil.virtual_memory().total / (1024 ** 3) * 0.9)
config = {'spark.driver.memory': f'{memory}g'}

hl.init(spark_conf=config, master=f'local[{args.cpus}]')

#Def from tpoterba on https://discuss.hail.is/t/importing-many-sample-specific-vcfs/2002
def multi_way_union_mts(mts: list, chunk_size: int) -> hl.MatrixTable:
    """Joins MatrixTables in the provided list
    :param list mts: list of MatrixTables to join together
    :param int chunk_size: number of MatrixTables to join per chunk
    :return: joined MatrixTable
    :rtype: MatrixTable
    """
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        print(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []
        for i in range(n_jobs):
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            print(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.missing(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )
            print(merged.aggregate((hl.agg.stats(hl.len(merged.__entries)), hl.len(merged.__cols))))
            next_stage.append(
                merged.checkpoint(
                    f"stage_{stage}_job_{i}.ht", overwrite=True
                )
            )
        print(f"done stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )

mts = [hl.read_matrix_table(table) for table in args.mt_files.split()]

mt_result = multi_way_union_mts(
    mts, 
    chunk_size = 64
)

mt_result.write(f'{args.output}/{args.species}.mt')
