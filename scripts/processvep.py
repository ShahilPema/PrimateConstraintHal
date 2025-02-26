import hail as hl
import psutil
import argparse
import os
import shutil

cpus = psutil.cpu_count(logical=True)

memory = int(psutil.virtual_memory().total/(1024 ** 3)*0.85)

config = {
    'spark.driver.memory': f'{memory}g',  #Set to total memory
    'spark.executor.memory': f'{memory}g'
}

hl.init(spark_conf=config, master=f'local[{cpus}]')

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process a VCF file with Hail.")
parser.add_argument("--vep_ht", required=True, help="Path to vep.ht.")
parser.add_argument("--output", required=True, help="Output directory")
args = parser.parse_args()

vep = hl.read_table(args.vep_ht)

vep = vep.repartition(cpus * 8)

vep = vep.checkpoint('A.ht', overwrite=True)

vep = hl.read_table('A.ht')

vep = vep.annotate(
    regions = vep.region.split(',')
)

vep = vep.explode('regions')

vep = vep.annotate(
    transcript = vep.regions.split('_')[0],
    region_type = vep.regions.split('_')[1],
    hg38_chr = vep.locus.contig,
    hg38_position = vep.locus.position
).drop('regions')

vep = vep.annotate(
    filtered_consequences = hl.if_else(
        vep.region_type == 'CDS',
        vep.vep.transcript_consequences.filter(
           lambda x: x.transcript_id == vep.transcript
        ),
        hl.missing('array<struct{allele_num: int32, cdna_end: int32, cdna_start: int32, gene_id: str, protein_end: int32, protein_start: int32, variant_allele: str, codons: str, consequence_terms: array<str>, strand: int32, transcript_id: str, impact: str, cds_start: int32, cds_end: int32, amino_acids: str, flags: array<str>, distance: int32, lof: str, lof_flags: str, lof_filter: str, lof_info: str, SpliceRegion: array<str>}>')
    )
)
#NOTE:SpliceRegion: array<str> may need to be replaced with spliceregion: array<str>, tssdistance: int32
vep = vep.annotate(
    codons = vep.filtered_consequences.codons,
    strand = vep.filtered_consequences.strand.first(),
    protein_start = vep.filtered_consequences.protein_start,
    cdna_start = vep.filtered_consequences.cdna_start,
    amino_acid = vep.filtered_consequences.amino_acids.first()
)

#{0: 2904693, 1: 404697456, None: 60666}
vep = vep.select('hg38_chr', 'hg38_position', 'codons', 'transcript', 'region_type', 'strand', 'protein_start', 'amino_acid')

vep = vep.annotate(
    codons = hl.if_else(
        (hl.is_defined(vep.codons) & (hl.len(vep.codons) > 0)),
        vep.codons[0],
        hl.missing(hl.tstr)
    )
)

vep = vep.annotate(
    codons = hl.if_else(
        hl.is_defined(vep.codons),
        vep.codons.strip().split("/")[0],
        hl.missing(hl.tstr)
    ),
    amino_acids = hl.if_else(
        hl.is_defined(vep.amino_acid),
        vep.amino_acid.strip().split("/")[0],
        hl.missing(hl.tstr)
    )
)

capital_bases = hl.set(['A', 'T', 'C', 'G'])

# Find index of the first match from our specific set
expr = hl.find(lambda x: capital_bases.contains(x[1]),
               hl.enumerate(vep.codons.split('')))[0]

vep = vep.annotate(
    codon_position=hl.if_else(
        hl.is_defined(vep.codons),
        expr,
        hl.missing(hl.tint32)
    ),
)

vep = vep.key_by()

vep = vep.select(
    'hg38_chr',
    'hg38_position',
    'codon_position',
    'codons',
    'transcript',
    'region_type',
    'strand',
    'protein_start',
    'amino_acid'
)

vep = vep.annotate(
    protein_start = hl.if_else(
        hl.is_defined(vep.protein_start) & (hl.len(vep.protein_start) > 0),
        vep.protein_start[0],
        hl.missing(hl.tint)
    )
)

vep = vep.checkpoint('B.ht', overwrite=True)

vep = vep.key_by('protein_start', 'transcript')

vep = vep.checkpoint('C.ht', overwrite=True)

vep_contiguous = vep.filter(vep.region_type == 'CDS')

vep_contiguous = vep_contiguous.group_by(
    'protein_start', 'transcript'
).aggregate(
    max_pos = hl.agg.max(vep_contiguous.hg38_position),
    min_pos = hl.agg.min(vep_contiguous.hg38_position)
)

vep_contiguous = vep_contiguous.annotate(
    contiguous = (vep_contiguous.max_pos - vep_contiguous.min_pos) < 3
    ).drop('max_pos', 'min_pos')

vep_contiguous = vep_contiguous.distinct()

vep_contiguous = vep_contiguous.checkpoint('D.ht', overwrite=True)

vep = vep.annotate(contiguous = vep_contiguous[vep.key].contiguous)

vep = vep.checkpoint('E.ht', overwrite=True)

vep = vep.key_by('hg38_chr', 'hg38_position', 'transcript', 'region_type')

vep = vep.distinct()

vep.write(f'{args.output}/vep_cdsinfo.ht', overwrite=True)

for file in ['A.ht', 'B.ht', 'C.ht', 'D.ht', 'E.ht']:
    if os.path.exists(file):
        shutil.rmtree(file)

##NOTE
##A few thousand positions are missing in the vep_annotated parquet because the mane select transcript was not included in the output. This behavior is replicated
##in the online VEP version. 
