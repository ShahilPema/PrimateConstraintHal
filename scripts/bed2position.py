import hail as hl
import psutil
import argparse
import shutil
import os
import gc
import math
import glob

parser = argparse.ArgumentParser(description="Perform liftover from species-specific BED to human BED using Hail.")
parser.add_argument("--species", required=True, help="Name of the species to process (e.g., Pithecia_pithecia).")
parser.add_argument("--human_fasta_path", required=True, help="Path to the human reference FASTA file.")
parser.add_argument("--human_index_path", required=True, help="Path to the index file for the human reference FASTA.")
parser.add_argument("--species_fasta_path", required=True, help="Path to the species-specific reference FASTA file.")
parser.add_argument("--species_index_path", required=True, help="Path to the index file for the species-specific reference FASTA.")
parser.add_argument("--species_bed_path", required=True, help="Path to the input species-specific BED file.")
parser.add_argument("--regionsxvar_table", required=True, help="Path to the regionsxvar Hail table (regionsxvar.ht).")
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

regionsxvar = hl.read_table(args.regionsxvar_table)

files = args.species_bed_path.split(' ')

lifttable = hl.import_table(files, no_header=True)

partitions = max(math.floor(len(files)/25), int(args.cpus)*2) #increase file size to 128MB but make sure each cpu has a partition
lifttable = lifttable.naive_coalesce(partitions)

lifttable = lifttable.annotate(
   human_coords = lifttable.f0,  # Original human coordinates
   species_coords = lifttable.f1,  # Target species coordinates  
   recip_human_coords = lifttable.f2   # Reciprocal human coordinates
)

# Create dictionary of species coords → reciprocal human coords
species_to_human = lifttable.group_by('species_coords').aggregate(
    recip_human_mappings = hl.agg.collect_as_set(lifttable.recip_human_coords)
).key_by('species_coords')

species_to_human = species_to_human.checkpoint('cpt1_1.ht', overwrite=True)

lifttable = lifttable.group_by('human_coords').aggregate(
    species_mappings = hl.agg.collect_as_set(lifttable.species_coords)
)

lifttable = lifttable.key_by()

lifttable = lifttable.checkpoint('cpt1.ht', overwrite=True)

human_rg = hl.ReferenceGenome.from_fasta_file(
    name="Homo_sapiens",
    fasta_file=args.human_fasta_path,
    index_file=args.human_index_path
)

species_rg = hl.ReferenceGenome.from_fasta_file(
    name=args.species,
    fasta_file=args.species_fasta_path,
    index_file=args.species_index_path
)

lifttable = lifttable.annotate(
    hg38_chr = lifttable.human_coords.split(':')[0],
    hg38_position = hl.int32(lifttable.human_coords.split(':')[1]),
    locus = hl.parse_locus(lifttable.human_coords, reference_genome=human_rg),
    human2species_mappings = hl.len(lifttable.species_mappings)
).drop('human_coords')

lifttable = lifttable.explode('species_mappings')

overlap = regionsxvar.group_by('hg38_chr', 'hg38_position').aggregate(
    overlapping_annotations = hl.int32(hl.agg.count())
)

overlap = overlap.checkpoint('cpt1_2.ht', overwrite=True)

lifttable = lifttable.annotate(
    species_contig = lifttable.species_mappings.split(':')[0],
    species_position = hl.int32(lifttable.species_mappings.split(':')[1]),
    species_locus = hl.parse_locus(lifttable.species_mappings, reference_genome=species_rg),
    reciprocal_mapping = species_to_human[lifttable.species_mappings].recip_human_mappings,
    overlapping_annotations = overlap[lifttable.hg38_chr, lifttable.hg38_position].overlapping_annotations
).drop('species_mappings')

lifttable = lifttable.checkpoint('cpt2_1.ht', overwrite=True)

#THIS STEP REQUIRES SOMEWHERE AROUND 10-20GB/CPU
lifttable = lifttable.annotate(
   hg38_ref_upstream = hl.get_sequence(lifttable.hg38_chr, lifttable.hg38_position, before = 15, reference_genome=human_rg).upper(),
   hg38_ref_downstream = hl.get_sequence(lifttable.hg38_chr, lifttable.hg38_position, after = 15, reference_genome=human_rg).upper(),
   hg38_pentamer = hl.get_sequence(lifttable.hg38_chr, lifttable.hg38_position,  before = 2, after = 2, reference_genome=human_rg).upper(),
   species_ref_upstream = hl.get_sequence(lifttable.species_contig, lifttable.species_position, before = 15, reference_genome=species_rg).upper(),
   species_ref_downstream = hl.get_sequence(lifttable.species_contig, lifttable.species_position, after = 15, reference_genome=species_rg).upper(),
)

lifttable = lifttable.checkpoint('cpt2_3.ht', overwrite=True)


lifttable = lifttable.annotate(
   species_ref_upstream_rev_comp = hl.reverse_complement(lifttable.species_ref_upstream),
   species_ref_downstream_rev_comp = hl.reverse_complement(lifttable.species_ref_downstream)
)

#99.9% of sequences are captured best by these comparisons when reverse and complement sequences are also considered
lifttable = lifttable.annotate(
   hammings = [
        hl.hamming(lifttable.hg38_ref_upstream, lifttable.species_ref_upstream),
        hl.hamming(lifttable.hg38_ref_downstream, lifttable.species_ref_downstream),
        hl.hamming(lifttable.hg38_ref_upstream, lifttable.species_ref_downstream_rev_comp),
        hl.hamming(lifttable.hg38_ref_downstream, lifttable.species_ref_upstream_rev_comp)
    ]
)

lifttable = lifttable.annotate(
   min_hamming = hl.min(lifttable.hammings)
)

comparisons = hl.array(["uO", "dO", "uO->dRC", "dO->uRC"])

lifttable = lifttable.checkpoint('cpt2_4.ht', overwrite=True)


lifttable = lifttable.annotate(
   min_hamming_orientation=hl.enumerate(lifttable.hammings)
       .filter(lambda idx_ham: idx_ham[1] == lifttable.min_hamming)  # Filter for min hamming values
       .map(lambda idx_ham: comparisons[idx_ham[0]])  # Map indexes to corresponding comparison labels
)  

plus_strand = hl.set([['dO'],['uO'],['uO', 'dO']])
minus_strand = hl.set([['dO->uRC'],['uO->dRC'],['uO->dRC', 'dO->uRC']])

lifttable = lifttable.annotate(
   strand_v_human_ref=hl.if_else(plus_strand.contains(lifttable.min_hamming_orientation),
                             hl.literal(1),
                             hl.if_else(minus_strand.contains(lifttable.min_hamming_orientation),
                                        hl.literal(-1),
                                        hl.missing(hl.tint)
                                       )
                            )
)

complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}

lifttable = lifttable.annotate(
   species_pentamer = hl.case() \
       .when(lifttable.strand_v_human_ref == 1, hl.get_sequence(lifttable.species_contig, lifttable.species_position,  before = 2, after = 2, reference_genome=species_rg).upper()) \
       .when(lifttable.strand_v_human_ref == -1, hl.reverse_complement(hl.get_sequence(lifttable.species_contig, lifttable.species_position,  before = 2, after = 2, reference_genome=species_rg).upper())) \
       .default(hl.missing(hl.tstr))
)

lifttable = lifttable.annotate(
   hg38_ref = lifttable.hg38_pentamer[2],
   species_refxstrand = lifttable.species_pentamer[2], 
   pentamer_error = hl.if_else(
       hl.is_defined(lifttable.hg38_pentamer),
       hl.hamming(lifttable.hg38_pentamer, lifttable.species_pentamer),
       hl.missing(hl.tint32)
   )
)

lifttable = lifttable.annotate(
   ref_match = hl.if_else(hl.is_defined(lifttable.species_refxstrand), 
                          hl.if_else(
                              lifttable.species_refxstrand == lifttable.hg38_ref,
                              hl.literal(1), 
                              hl.literal(0)),
                          hl.missing(hl.tint32)
                         ),
   species_ref = hl.case() \
       .when(lifttable.strand_v_human_ref == 1, lifttable.species_refxstrand) \
       .when(lifttable.strand_v_human_ref == -1, lifttable.species_refxstrand.translate(complement)) \
       .default(hl.missing(hl.tstr))
)

lifttable = lifttable.checkpoint('cpt2.ht', overwrite=True)

pos_table = lifttable.annotate(
       sequence_info = hl.struct(
           ref_match = lifttable.ref_match,
           strand_v_human_ref = lifttable.strand_v_human_ref,
           min_hamming = hl.if_else(
               hl.is_defined(lifttable.strand_v_human_ref),
               lifttable.min_hamming,
               hl.missing(hl.tint32)
           ),
           species_pentamer = lifttable.species_pentamer,
           hg38_pentamer = lifttable.hg38_pentamer,
           pentamer_error = lifttable.pentamer_error,
           reciprocal_mapping = lifttable.reciprocal_mapping,
           human2species_mappings = lifttable.human2species_mappings,
           overlapping_annotations = lifttable.overlapping_annotations
       )
).rename(
   {'sequence_info': f'{args.species}_pos_info', 'species_locus':f'locus_{args.species}'}
)

pos_table = pos_table.checkpoint('cpt3.ht', overwrite=True)

pos_table = pos_table.key_by(pos_table[f'locus_{args.species}'])

pos_table = pos_table.select(
   pos_table.locus,
   f'{args.species}_pos_info'
)

pos_table.write(f'{args.output}/{args.species}_pos_info.ht', overwrite=True)