import hail as hl
import psutil
import argparse
import shutil
import os
import gc
import math

parser = argparse.ArgumentParser(description="Perform liftover from species-specific BED to human BED using Hail.")
parser.add_argument("--species", required=True, help="Name of the species to process (e.g., Pithecia_pithecia).")
parser.add_argument("--sample_mts", required=True, help="Sample hail tables to be merged ")
parser.add_argument("--full_ht", required=True, help="Species variant info hail table.")
parser.add_argument("--cpus", required=True, help="cpus for task")
parser.add_argument("--memory", required=True, help="memory for task")
parser.add_argument("--contig_count", required=True, help="Total contigs.")
parser.add_argument("--tmpdir", required=True, help="Temporary directory")
parser.add_argument("--outdir", required=True, help="Total system cpus")
args = parser.parse_args()

#NOTE change cpus so each partition is more than 128MB. 

config = {
    'spark.driver.memory': f'{args.memory}g',  #Set to total memory
    'spark.executor.memory': f'{args.memory}g',
    'spark.local.dir': args.tmpdir
}

hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tempdir=args.tmpdir)

full_ht = hl.read_table(args.full_ht)

full_ht = full_ht.annotate(
    contig = full_ht[f'{args.species}_locus'].contig,
    position = full_ht[f'{args.species}_locus'].position
)

contig_count = int(args.contig_count)
if contig_count > 100000: 
	full_ht = full_ht.key_by('contig', 'position', f'{args.species}_alleles')

sample_mt_list = args.sample_mts.split()
column_data = None
count = 1

for file in sample_mt_list:
    print(file)
    mt = hl.read_matrix_table(file)
    mt = mt.key_cols_by()
    if column_data == None:
        column_data = mt.cols()
    else:
        column_data=column_data.union(mt.cols())
    sample_id = mt.s.collect()[0]
    mt = mt.select_cols()
    ht = mt.entries()
    if contig_count > 100000:
        ht = ht.annotate(
            contig = ht.locus.contig,
            position = ht.locus.position
        )
        ht = ht.rename({'alleles': f'{args.species}_alleles'})
        ht = ht.key_by('contig', 'position', f'{args.species}_alleles')
    else:
        ht = ht.rename({'locus': f'{args.species}_locus', 'alleles': f'{args.species}_alleles'})
        ht = ht.key_by(f'{args.species}_locus', f'{args.species}_alleles')
    full_ht = full_ht.annotate(GT = ht[full_ht.key].GT)
    full_ht = full_ht.annotate(GT = hl.case() \
                                .when(hl.is_missing(full_ht.alleles), hl.missing(hl.tcall)) \
                                .when(hl.is_defined(full_ht.GT), full_ht.GT) \
                                .default(hl.call(0,0))
                          )
    full_ht = full_ht.rename({'GT':sample_id})
    full_ht = full_ht.checkpoint(f'./tmp/checkpoint_{count}.ht', overwrite=True)
    count += 1

full_ht = full_ht.annotate(
    species_data = hl.struct(
        species_locus = full_ht[f'{args.species}_locus'],
        species_alleles = full_ht[f'{args.species}_alleles'],
        transcript_mappings = full_ht[f'{args.species}_transcript_mapping']
    )
)

full_ht = full_ht.rename({'species_data': f'{args.species}_data'})

full_ht = full_ht.checkpoint('./tmp/A.ht', overwrite=True)

full_ht = full_ht.drop(f'{args.species}_locus', f'{args.species}_alleles', f'{args.species}_transcript_mapping', 'contig', 'position')

full_ht.write(f'{args.outdir}/{args.species}_lifted_vars.ht', overwrite=True)

column_data.write(f'{args.outdir}/{args.species}_colinfo_vars.ht', overwrite=True)

if os.path.isdir('tmp'):
    shutil.rmtree('tmp')