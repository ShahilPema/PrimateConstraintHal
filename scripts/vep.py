import hail as hl
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process a VCF file with Hail.")
parser.add_argument("--cpus", required=True, help="Total cpus on machine")
parser.add_argument("--memory", required=True, help="Total memory on machine")
parser.add_argument("--regionsxvar_table", required=True, help="Path to table of positions to be annotated with VEP.")
parser.add_argument("--vep_config", required=True, help="Config file for VEP.")
parser.add_argument("--tmpdir", required=True, help="Temporary directory")
parser.add_argument("--output", required=True, help="Output directory")
args = parser.parse_args()

config = {
    'spark.driver.memory': f'{args.memory}',  #Set to total memory
    'spark.executor.memory': f'{args.memory}',
    'spark.local.dir': args.tmpdir,
    'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}',
    'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}'
}

hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)

bed_ht = hl.read_table(args.regionsxvar_table)
vep_annotate = bed_ht.filter(bed_ht.regions.contains('-CDS'))
vep_annotate = vep_annotate.checkpoint('cpt1.ht', overwrite=True)

if vep_annotate.count() > 0:
    vep_annotate = hl.vep(vep_annotate, args.vep_config)
    vep_annotate = vep_annotate.checkpoint('cpt2.ht', overwrite=True)
    bed_ht = bed_ht.annotate(
        vep = vep_annotate[bed_ht.key]
    )
    bed_ht.write(f'{args.output}/vep.ht', overwrite=True)
else:
    vep_col = 'struct{input: str, id: str, end: int32, transcript_consequences: array<struct{allele_num: int32, cdna_end: int32, cdna_start: int32, gene_id: str, protein_end: int32, protein_start: int32, variant_allele: str, codons: str, consequence_terms: array<str>, strand: int32, transcript_id: str, impact: str, cds_start: int32, cds_end: int32, amino_acids: str, flags: array<str>, distance: int32, lof: str, lof_flags: str, lof_filter: str, lof_info: str, spliceregion: array<str>, tssdistance: int32}>, intergenic_consequences: array<struct{allele_num: int32, variant_allele: str, consequence_terms: array<str>, impact: str}>, motif_feature_consequences: array<struct{allele_num: int32, cdna_end: int32, cdna_start: int32, gene_id: str, protein_end: int32, protein_start: int32, variant_allele: str, codons: str, consequence_terms: array<str>, strand: int32, transcript_id: str, impact: str, cds_start: int32, cds_end: int32, amino_acids: str, flags: array<str>}>, regulatory_feature_consequences: array<struct{allele_num: int32, cdna_end: int32, cdna_start: int32, gene_id: str, protein_end: int32, protein_start: int32, variant_allele: str, codons: str, consequence_terms: array<str>, strand: int32, transcript_id: str, impact: str, cds_start: int32, cds_end: int32, amino_acids: str, flags: array<str>, biotype: str, regulatory_feature_id: str}>, assembly_name: str, seq_region_name: str, most_severe_consequence: str, start: int32, allele_string: str, strand: int32}'
    vep_annotate = vep_annotate.annotate(
        vep = hl.missing(vep_col)
    )
    vep_annotate.write(f'{args.output}/vep.ht', overwrite=True)

##NOTE
##Some positions are missing in the vep_annotated parquet because the mane select transcript was not included in the output. This behavior is replicated
##in the online VEP version.
# vep = pl.scan_parquet('vep.parquet/part-00000-5b571740-bf37-4392-8e7c-adf58df8cd93-c000.snappy.parquet').select(
#     'hg38_chr','hg38_position').with_columns(pl.lit(1).alias('vep_annotated'))

# missing = pl.scan_parquet('/home/ubuntu/PrimateConstraintHal/127vars.parquet').with_columns( pl.col('hg38_position').cast(pl.Int32).alias('hg38_position'))

# missing = missing.join(vep, how='left', on=['hg38_chr', 'hg38_position'], coalesce=True).collect()

# missing = missing.filter(pl.col('vep_annotated').is_null())
