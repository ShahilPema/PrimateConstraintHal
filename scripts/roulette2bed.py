import hail as hl
import psutil
import argparse
import shutil
import os
import gc
import math
import glob

parser = argparse.ArgumentParser(description="Perform liftover from species-specific BED to human BED using Hail.")
parser.add_argument("--roulette_table", required=True, help="Path to the roulette Hail table (roulette.ht).")
parser.add_argument("--cpus", required=True, help="Cpus for task")
parser.add_argument("--memory", required=True, help="Memory for task")
parser.add_argument("--output", required=True, help="Directory where the output files will be written.")
parser.add_argument("--tmpdir", required=True, help="Location for tmp file storage")
args = parser.parse_args()

#NOTE change cpus so each partition is more than 128MB. 
#The resulting tables now include all positions that a given human position maps to in the target species.

config = {
    'spark.driver.memory': f'{args.memory}',  #Set to total memory
    'spark.executor.memory': f'{args.memory}',
    'spark.local.dir': args.tmpdir,
    'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}',
    'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}'
}

hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)

roulette = hl.read_table(args.roulette_table)

roulette = roulette.filter(roulette.roulette_data.roulette_nc_scaling_site == True)

roulette = roulette.key_by('locus')

roulette = roulette.select()

roulette = roulette.distinct()

roulette = roulette.key_by()

roulette = roulette.annotate(
    contig = hl.str(roulette.locus.contig),
    start = hl.int(roulette.locus.position) - 1,
    end = hl.int(roulette.locus.position),
)

roulette = roulette.drop('locus')

roulette.export(args.output, header=False)