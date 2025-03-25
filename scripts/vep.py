import hail as hl
from pyspark.sql import SparkSession
import os
import argparse
import psutil

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process a VCF file with Hail.")
parser.add_argument("--cpus", required=True, help="Total cpus on machine")
parser.add_argument("--position_table", required=True, help="Path to table of positions to be annotated with VEP.")
parser.add_argument("--human_fasta_path", required=True, help="Path to the FASTA file for the primate genome.")
parser.add_argument("--human_index_path", required=True, help="Path to the FASTA index file (.fai).")
parser.add_argument("--vep_config", required=True, help="Config file for VEP.")
parser.add_argument("--output", required=True, help="Output directory")
args = parser.parse_args()

memory = int(psutil.virtual_memory().total / (1024 ** 3) * 0.9)
config = {'spark.driver.memory': f'{memory}g'}

hl.init(spark_conf=config, master=f'local[{args.cpus}]')

rg = hl.ReferenceGenome.from_fasta_file(
    name='Homo_sapiens',
    fasta_file=args.human_fasta_path,
    index_file=args.human_index_path
)

spark = SparkSession.builder.appName("ParquetReader").getOrCreate()

bed_data = spark.read.parquet(args.position_table)
bed_ht = hl.Table.from_spark(bed_data)
bed_ht = bed_ht.repartition(int(args.cpus))
bed_ht = bed_ht.annotate(variant_str = bed_ht.hg38_chr + ":" + hl.str(bed_ht.hg38_position))
bed_ht = bed_ht.annotate(locus = hl.parse_locus(bed_ht.variant_str,reference_genome='Homo_sapiens'))
bed_ht = bed_ht.annotate(
    ref=bed_ht.locus.sequence_context()
)
bed_ht = bed_ht.annotate(
    possible_alts=hl.if_else(
        bed_ht.ref != "N",
        hl.set(['A', 'T', 'C', 'G']).remove(bed_ht.ref),
        hl.set(['N'])
    )
)
bed_ht = bed_ht.explode('possible_alts')
bed_ht = bed_ht.annotate(
    alleles=hl.array([bed_ht.ref, bed_ht.possible_alts])
)

bed_ht = bed_ht.key_by(bed_ht.locus, bed_ht.alleles)

bed_ht = bed_ht.select(
    hg38_chr = bed_ht.hg38_chr,
    hg38_position = hl.int32(bed_ht.hg38_position),
    REF = bed_ht.ref,
    ALT = bed_ht.possible_alts,
    region = bed_ht.region
)

bed_ht = bed_ht.checkpoint(f'{args.output}/vep_cpt.ht', overwrite=True)
bed_ht = hl.vep(bed_ht, args.vep_config)
bed_ht.write(f'{args.output}/vep.ht', overwrite=True)

##NOTE
##Some positions are missing in the vep_annotated parquet because the mane select transcript was not included in the output. This behavior is replicated
##in the online VEP version.
# vep = pl.scan_parquet('vep.parquet/part-00000-5b571740-bf37-4392-8e7c-adf58df8cd93-c000.snappy.parquet').select(
#     'hg38_chr','hg38_position').with_columns(pl.lit(1).alias('vep_annotated'))

# missing = pl.scan_parquet('/home/ubuntu/PrimateConstraintHal/127vars.parquet').with_columns( pl.col('hg38_position').cast(pl.Int32).alias('hg38_position'))

# missing = missing.join(vep, how='left', on=['hg38_chr', 'hg38_position'], coalesce=True).collect()

# missing = missing.filter(pl.col('vep_annotated').is_null())
