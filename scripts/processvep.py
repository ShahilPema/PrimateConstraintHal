import hail as hl
import psutil
import argparse
import os
import shutil

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process a VCF file with Hail.")
parser.add_argument("--vep_ht", required=True, help="Path to vep.ht.")
parser.add_argument("--memory", required=True, help="Memory available for this task.")
parser.add_argument("--cpus", required=True, help="Cpus available for this task.")
parser.add_argument("--output", required=True, help="Output directory")
parser.add_argument("--tmpdir", required=True, help="Local temp dir")
args = parser.parse_args()

config = {
    'spark.driver.memory': f'{args.memory}',  #Set to total memory
    'spark.executor.memory': f'{args.memory}',
    'spark.local.dir': args.tmpdir,
    'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}',
    'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}'
}

hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)

vep = hl.read_table(args.vep_ht)

vep = vep.repartition(int(args.cpus) * 4)

vep = vep.checkpoint('A.ht', overwrite=True)

vep = vep.annotate(
    regions = vep.regions.split(',')
)

vep = vep.explode('regions')

vep = vep.annotate(
    transcript = vep.regions.split('-')[0],
    region_type = vep.regions.split('-')[1],
    hg38_chr = vep.locus.contig,
    hg38_position = vep.locus.position
).drop('regions')

vep = vep.annotate(
    filtered_consequences = hl.if_else(
        vep.region_type == 'CDS',
        vep.vep.transcript_consequences.filter(
           lambda x: x.transcript_id == vep.transcript
        ),
        hl.missing(vep.vep.transcript_consequences.dtype)
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

"""
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
"""

vep = vep.checkpoint('E.ht', overwrite=True)

vep = vep.key_by('hg38_chr', 'hg38_position', 'transcript', 'region_type')

vep = vep.distinct()

vep.write(f'{args.output}/vep_cdsinfo.ht', overwrite=True)

for file in ['A.ht', 'B.ht', 'C.ht', 'D.ht', 'E.ht']:
    if os.path.exists(file):
        shutil.rmtree(file)
