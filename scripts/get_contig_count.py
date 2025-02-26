#!/usr/bin/env python3
import argparse
import hail as hl

parser = argparse.ArgumentParser(description="Get number of unique contigs in genome.")
parser.add_argument("--species", required=True, help="Species.")
parser.add_argument("--species_fasta_path", required=True, help="Path to the species-specific reference FASTA file.")
parser.add_argument("--species_index_path", required=True, help="Path to the index file for the species-specific reference FASTA.")
args = parser.parse_args()

hl.init()

species_rg = hl.ReferenceGenome.from_fasta_file(
    name=args.species,
    fasta_file=args.species_fasta_path, 
    index_file=args.species_index_path
)

print(len(species_rg.contigs))
