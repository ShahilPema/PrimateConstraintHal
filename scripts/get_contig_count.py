#!/usr/bin/env python3
import argparse
import hail as hl

parser = argparse.ArgumentParser(description="Get number of unique contigs in genome.")
parser.add_argument("--cpus", required=True, help="Total cpus on machine")
parser.add_argument("--memory", required=True, help="Total memory on machine")
parser.add_argument("--species", required=True, help="Species.")
parser.add_argument("--species_fasta_path", required=True, help="Path to the species-specific reference FASTA file.")
parser.add_argument("--species_index_path", required=True, help="Path to the index file for the species-specific reference FASTA.")
parser.add_argument("--tmpdir", required=True, help="Temporary directory")
args = parser.parse_args()

config = {
    'spark.driver.memory': f'{args.memory}',  #Set to total memory
    'spark.executor.memory': f'{args.memory}',
    'spark.local.dir': args.tmpdir,
    'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}',
    'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}'
}
 
hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)

species_rg = hl.ReferenceGenome.from_fasta_file(
    name=args.species,
    fasta_file=args.species_fasta_path, 
    index_file=args.species_index_path
)

print(len(species_rg.contigs))
