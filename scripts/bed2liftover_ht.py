import hail as hl
import psutil
import argparse
import shutil
import os
import gc
import math

parser = argparse.ArgumentParser(description="Perform liftover from species-specific BED to human BED using Hail.")
parser.add_argument("--species", required=True, help="Name of the species to process (e.g., Pithecia_pithecia).")
parser.add_argument("--human_fasta_path", required=True, help="Path to the human reference FASTA file.")
parser.add_argument("--human_index_path", required=True, help="Path to the index file for the human reference FASTA.")
parser.add_argument("--species_fasta_path", required=True, help="Path to the species-specific reference FASTA file.")
parser.add_argument("--species_index_path", required=True, help="Path to the index file for the species-specific reference FASTA.")
parser.add_argument("--species_bed_path", required=True, help="Path to the input species-specific BED file.")
parser.add_argument("--vep_annotations", required=True, help="Path to the VEP annotations Hail table (vep.ht).")
parser.add_argument("--all_cpus", required=True, help="Total system cpus")
parser.add_argument("--output", required=True, help="Directory where the output files will be written.")
args = parser.parse_args()

#NOTE change cpus so each partition is more than 128MB. 

hl_cpus = 90

forks = math.floor(int(args.all_cpus)/hl_cpus)

memory = int(psutil.virtual_memory().total/(1024 ** 3)*0.85/forks)

config = {
    'spark.driver.memory': f'{memory}g',  #Set to total memory
    'spark.executor.memory': f'{memory}g'
}

hl.init(spark_conf=config, master=f'local[{hl_cpus}]')

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

lifttable = hl.import_table(args.species_bed_path, no_header=True, min_partitions = hl_cpus * 2)

lifttable = lifttable.annotate(
   hg38_chr = lifttable.f0,
   hg38_position = hl.int32(lifttable.f2),
   species_contig = lifttable.f3.split(':')[0],
   species_position = hl.int32(lifttable.f3.split(':')[1])
)
lifttable = lifttable.select('hg38_chr', 'hg38_position', 'species_contig', 'species_position')

lifttable = lifttable.annotate(
    locus = hl.parse_locus(lifttable.hg38_chr + ":" + hl.str(lifttable.hg38_position), reference_genome="Homo_sapiens"),
    species_locus = hl.parse_locus(lifttable.species_contig + ":" + hl.str(lifttable.species_position), reference_genome=args.species)
)

lifttable = lifttable.annotate(
   species_ref_upstream = hl.get_sequence(lifttable.species_contig, lifttable.species_position, before = 15, reference_genome=args.species).upper(),
   hg38_ref_upstream = hl.get_sequence(lifttable.hg38_chr, lifttable.hg38_position, before = 15, reference_genome="Homo_sapiens").upper(),
   species_ref_downstream = hl.get_sequence(lifttable.species_contig, lifttable.species_position, after = 15, reference_genome=args.species).upper(),
   hg38_ref_downstream = hl.get_sequence(lifttable.hg38_chr, lifttable.hg38_position, after = 15, reference_genome="Homo_sapiens").upper(),
   hg38_pentamer = hl.get_sequence(lifttable.hg38_chr, lifttable.hg38_position,  before = 2, after = 2, reference_genome="Homo_sapiens").upper(),
)

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

lifttable = lifttable.annotate(
   min_hamming_orientation=hl.enumerate(lifttable.hammings)
       .filter(lambda idx_ham: idx_ham[1] == lifttable.min_hamming)  # Filter for min hamming values
       .map(lambda idx_ham: comparisons[idx_ham[0]])  # Map indexes to corresponding comparison labels
)  

plus_strand = hl.set([['dO'],['uO'],['uO', 'dO']])
minus_strand = hl.set([['dO->uRC'],['uO->dRC'],['uO->dRC', 'dO->uRC']])

lifttable = lifttable.annotate(
   strand_v_human=hl.if_else(plus_strand.contains(lifttable.min_hamming_orientation),
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
       .when(lifttable.strand_v_human == 1, hl.get_sequence(lifttable.species_contig, lifttable.species_position,  before = 2, after = 2, reference_genome=args.species).upper()) \
       .when(lifttable.strand_v_human == -1, hl.reverse_complement(hl.get_sequence(lifttable.species_contig, lifttable.species_position,  before = 2, after = 2, reference_genome=args.species).upper())) \
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
       .when(lifttable.strand_v_human == 1, lifttable.species_refxstrand) \
       .when(lifttable.strand_v_human == -1, lifttable.species_refxstrand.translate(complement)) \
       .default(hl.missing(hl.tstr))
)

lifttable = lifttable.checkpoint('A.ht', overwrite=True)

pos_table = lifttable.annotate(
       sequence_info = hl.struct(
           ref_match = lifttable.ref_match,
           strand_v_human = lifttable.strand_v_human,
           min_hamming = lifttable.min_hamming,
           species_pentamer = lifttable.species_pentamer,
           hg38_pentamer = lifttable.hg38_pentamer,
           pentamer_error = lifttable.pentamer_error,
       )
).rename(
   {'sequence_info': f'{args.species}_pos_info', 'species_locus':f'locus_{args.species}'}
)
pos_table = pos_table.key_by(pos_table[f'locus_{args.species}'])

pos_table = pos_table.select(
   pos_table.locus,
   f'{args.species}_pos_info'
)

pos_table.write(f'{args.output}/{args.species}_pos_info.ht', overwrite=True)

del(pos_table)

gc.collect()

lifttable = lifttable.select('hg38_chr', 'hg38_position', 'species_contig', 'species_position', 'species_pentamer', 'strand_v_human',
                           'hg38_ref', 'species_ref', 'ref_match', 'locus', 'species_locus')

translate = hl.literal({
    'TTT': 'F', 'TTC': 'F',  # Phenylalanine
    'TTA': 'L', 'TTG': 'L',  # Leucine
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',  # Leucine
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',  # Isoleucine
    'ATG': 'M',  # Methionine (Start codon)
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # Valine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',  # Serine
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  # Proline
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  # Threonine
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  # Alanine
    'TAT': 'Y', 'TAC': 'Y',  # Tyrosine
    'TAA': '*', 'TAG': '*', 'TGA': '*',  # Stop codons
    'CAT': 'H', 'CAC': 'H',  # Histidine
    'CAA': 'Q', 'CAG': 'Q',  # Glutamine
    'AAT': 'N', 'AAC': 'N',  # Asparagine
    'AAA': 'K', 'AAG': 'K',  # Lysine
    'GAT': 'D', 'GAC': 'D',  # Aspartic acid
    'GAA': 'E', 'GAG': 'E',  # Glutamic acid
    'TGT': 'C', 'TGC': 'C',  # Cysteine
    'TGG': 'W',  # Tryptophan
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',  # Arginine
    'AGT': 'S', 'AGC': 'S',  # Serine
    'AGA': 'R', 'AGG': 'R',  # Arginine
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'   # Glycine
})


vep = hl.read_table(args.vep_annotations)
vep = vep.key_by('hg38_chr', 'hg38_position')

overlap = vep.group_by('hg38_chr', 'hg38_position').aggregate(
    overlap = hl.agg.count()
)

vep = vep.annotate(
    overlap=overlap[vep.key].overlap
)

lifttable = lifttable.key_by('hg38_chr', 'hg38_position').join(vep)

lifttable = lifttable.checkpoint('B1.ht', overwrite=True)

subset = lifttable.filter(lifttable.contiguous == False)

subset = subset.repartition(hl_cpus * 2)

subset = subset.group_by(subset.transcript, subset.region_type, subset.protein_start).aggregate(
    species_loci = hl.zip(hl.agg.collect(subset.species_contig), hl.agg.collect(subset.species_position))
)

subset = subset.annotate(
    sorted_species_loci = hl.sorted(subset.species_loci, key=lambda x: x[1])  # Sort by position (index 1)
)

subset = subset.annotate(
    species_ref_codon_noncontiguous = hl.str('').join(
        hl.map(
            lambda locus: hl.get_sequence(locus[0], locus[1], reference_genome=args.species).upper(),
            subset.sorted_species_loci
        )
    )
)

subset = subset.select('species_ref_codon_noncontiguous')

subset = subset.checkpoint('subset.ht', overwrite=True)

lifttable =lifttable.key_by('transcript', 'region_type', 'protein_start')

lifttable = lifttable.annotate(
    species_codon_strand = hl.if_else(
        hl.is_defined(lifttable.strand),
        lifttable.strand_v_human * lifttable.strand,
        hl.missing(hl.tint)
    )
)


lifttable = lifttable.annotate(
    species_codon_contiguous = hl.case() \
        .when((lifttable.strand == 1), lifttable.species_pentamer[2-lifttable.codon_position:5-lifttable.codon_position])
        .when((lifttable.strand == -1), hl.reverse_complement(lifttable.species_pentamer[lifttable.codon_position:lifttable.codon_position+3]))
        .default(hl.missing(hl.tstr)),
    species_codon_noncontiguous = hl.if_else(
        lifttable.species_codon_strand == 1,
        subset[lifttable.key].species_ref_codon_noncontiguous,
        hl.if_else(
            lifttable.species_codon_strand == -1,
            hl.reverse_complement(subset[lifttable.key].species_ref_codon_noncontiguous),
            hl.missing(hl.tstr)
            )
        ),
    human_codon = lifttable.codons.upper()
)

lifttable = lifttable.annotate(
    species_codon = hl.if_else(
        hl.is_defined(lifttable.species_codon_contiguous),
        lifttable.species_codon_contiguous,
        hl.if_else(
            hl.is_defined(lifttable.species_codon_noncontiguous),
            lifttable.species_codon_noncontiguous,
            hl.missing(hl.tstr)
        )
    )
)

lifttable = lifttable.annotate(
    codon_match = hl.if_else(
        hl.is_defined(lifttable.human_codon),
        hl.if_else(
            lifttable.human_codon == lifttable.species_codon,
            hl.literal(1),
            hl.literal(0),
            ),
        hl.missing(hl.tint)
    ),
    human_aa = hl.if_else(
        hl.is_defined(lifttable.human_codon),
        translate.get(lifttable.human_codon),
        hl.missing(hl.tstr)
    ),
    species_aa = hl.if_else(
        hl.is_defined(lifttable.species_codon),
        translate.get(lifttable.species_codon),
        hl.missing(hl.tstr)
    )
)

lifttable = lifttable.annotate(
    aa_match = hl.if_else(
        hl.is_defined(lifttable.human_codon),
        hl.if_else(
            lifttable.human_aa == lifttable.species_aa,
            hl.literal(1),
            hl.literal(0),
            ),
        hl.missing(hl.tint)
    )
)


lifttable = lifttable.annotate(
    possible_alts = hl.if_else(
        lifttable.hg38_ref != "N",
        hl.set(['A', 'T', 'C', 'G']).remove(lifttable.hg38_ref),
        hl.set(['N'])
    )
)

lifttable = lifttable.checkpoint('B.ht', overwrite=True)

lifttable = lifttable.explode('possible_alts')
lifttable = lifttable.annotate(
    alleles=hl.array([lifttable.hg38_ref, lifttable.possible_alts])
)

lifttable = lifttable.annotate(
    human_alt_codon = hl.if_else(
        hl.is_defined(lifttable.codon_position) & (lifttable.strand == 1),
        lifttable.human_codon[0:lifttable.codon_position] + lifttable.possible_alts + lifttable.human_codon[lifttable.codon_position + 1:],
        hl.if_else(
            hl.is_defined(lifttable.codon_position) & (lifttable.strand == -1),
            lifttable.human_codon[0:lifttable.codon_position] + lifttable.possible_alts.translate(complement) + lifttable.human_codon[lifttable.codon_position + 1:],
            hl.missing(hl.tstr)
        )
    ),
    species_alt_codon = hl.if_else(
        hl.is_defined(lifttable.codon_position) & (lifttable.strand == 1),
        lifttable.species_codon[0:lifttable.codon_position] + lifttable.possible_alts + lifttable.species_codon[lifttable.codon_position + 1:],
        hl.if_else(
            hl.is_defined(lifttable.codon_position) & (lifttable.strand == -1),
            lifttable.species_codon[0:lifttable.codon_position] + lifttable.possible_alts.translate(complement) + lifttable.species_codon[lifttable.codon_position + 1:],
            hl.missing(hl.tstr)
        )
    )
)

lifttable = lifttable.annotate(
    human_alt_aa = translate.get(lifttable.human_alt_codon),
    species_alt_aa = translate.get(lifttable.species_alt_codon)
)

def annotate_alt_consequence(ref_aa, alt_aa, codon_position, codon_match=None, ref_match=None, species=False):
    return hl.case() \
        .when((species == True) & (ref_match == 0), "ref_mismatch") \
        .when((ref_aa == "M") & (alt_aa != "M") & (codon_position == 1), "start_lost") \
        .when((ref_aa != alt_aa) & (ref_aa != "*") & (alt_aa != "*"), "missense_variant") \
        .when((ref_aa == "*") & (alt_aa == "*"), "stop_retained") \
        .when((ref_aa != "*") & (alt_aa == "*"), "stop_gained") \
        .when((ref_aa == "*") & (alt_aa != "*"), "stop_lost") \
        .when(ref_aa == alt_aa, "synonymous_variant") \
        .when(hl.is_missing(ref_aa), hl.missing(hl.tstr)) \
        .default("other")

# Apply the function for both human and species annotations
lifttable = lifttable.annotate(
    human_alt_cons=annotate_alt_consequence(lifttable.human_aa, lifttable.human_alt_aa, lifttable.codon_position),
    species_alt_cons=annotate_alt_consequence(lifttable.species_aa, lifttable.species_alt_aa, lifttable.codon_position, lifttable.codon_match, lifttable.ref_match, species=True),
    codon_mismatch_cons = hl.if_else(
        lifttable.codon_match == 0,
        annotate_alt_consequence(lifttable.human_aa, lifttable.species_aa, lifttable.codon_position),
        hl.missing(hl.tstr)
    )
)

lifttable = lifttable.annotate(
    alt_cons_match = hl.case() \
         .when((lifttable.species_alt_cons == 'ref_mismatch'), hl.missing(hl.tint))
         .when((lifttable.human_alt_cons == lifttable.species_alt_cons) & (hl.is_defined(lifttable.human_alt_cons)), hl.literal(1)) \
         .when((lifttable.human_alt_cons != lifttable.species_alt_cons) & (hl.is_defined(lifttable.human_alt_cons)), hl.literal(0)) \
         .default(hl.missing(hl.tint))
)

lifttable = lifttable.key_by()

def process_alleles(alleles: hl.tarray, species_ref: hl.tstr):
    ref_original = alleles[0]
    alt_original = alleles[1]
    return hl.case() \
        .when(ref_original == "N", hl.array([species_ref, species_ref])) \
        .when(alt_original == species_ref, hl.array([species_ref, ref_original])) \
        .default(hl.array([species_ref, alt_original]))

lifttable = lifttable.annotate(
    species_alleles = hl.case() \
        .when(lifttable.ref_match == 1,
              hl.case() \
              .when(lifttable.strand_v_human == 1,lifttable.alleles) \
              .when(lifttable.strand_v_human == -1, hl.map(lambda x: x.translate(complement), lifttable.alleles)) \
              .default(hl.missing(hl.tarray(hl.tstr)))
             )
        .when(lifttable.ref_match == 0,
              hl.case() \
              .when(lifttable.strand_v_human == 1, process_alleles(lifttable.alleles, lifttable.species_ref)) \
              .when(lifttable.strand_v_human == -1, hl.map(lambda x: x.translate(complement), process_alleles(lifttable.alleles, lifttable.species_ref))) \
              .default(hl.missing(hl.tarray(hl.tstr)))
             )
        .default(hl.missing(hl.tarray(hl.tstr))),
)

lifttable = lifttable.checkpoint('C.ht', overwrite=True)

lifttable_no_overlap = lifttable.filter(lifttable.overlap == 1)

lifttable_no_overlap = lifttable_no_overlap.annotate(
    transcript_mapping = hl.set([hl.struct(
            transcript = lifttable_no_overlap.transcript,
            region_type = lifttable_no_overlap.region_type,
            strand = lifttable_no_overlap.strand,
            codon_match = lifttable_no_overlap.codon_match,
            aa_match = lifttable_no_overlap.aa_match,
            codon_mismatch_cons = lifttable_no_overlap.codon_mismatch_cons,
            human_alt_cons = lifttable_no_overlap.human_alt_cons,
            species_alt_cons = lifttable_no_overlap.species_alt_cons,
            alt_cons_match = lifttable_no_overlap.alt_cons_match
        )])
)

lifttable_no_overlap = lifttable_no_overlap.select('locus', 'alleles', 'species_locus', 'species_alleles', 'transcript_mapping')

lifttable_overlap = lifttable.filter(lifttable.overlap > 1)
lifttable_overlap = lifttable_overlap.group_by('locus', 'alleles').aggregate(
    species_locus = hl.agg.collect(lifttable_overlap.species_locus)[0],
    species_alleles = hl.agg.collect(lifttable_overlap.species_alleles)[0],
    transcript_mapping=hl.agg.collect_as_set(
        hl.struct(
            transcript = lifttable_overlap.transcript,
            region_type = lifttable_overlap.region_type,
            strand = lifttable_overlap.strand,
            codon_match = lifttable_overlap.codon_match,
            aa_match = lifttable_overlap.aa_match,
            codon_mismatch_cons = lifttable_overlap.codon_mismatch_cons,
            human_alt_cons = lifttable_overlap.human_alt_cons,
            species_alt_cons = lifttable_overlap.species_alt_cons,
            alt_cons_match = lifttable_overlap.alt_cons_match
        )
    )
)

lifttable_overlap = lifttable_overlap.key_by()


lifttable = lifttable_overlap.union(lifttable_no_overlap)

lifttable = lifttable.annotate(
    alleles = hl.if_else(
        hl.is_defined(lifttable.species_alleles),
        lifttable.alleles,
        hl.missing(hl.tarray(hl.tstr))
    )
)

lifttable = lifttable.rename({ 'species_locus' :f'{args.species}_locus', 'species_alleles' : f'{args.species}_alleles', 'transcript_mapping' : f'{args.species}_transcript_mapping'})
lifttable.write(f'{args.output}/{args.species}_var_info.ht', overwrite=True)

for file in ['A.ht', 'B.ht', 'B1.ht', 'subset.ht', 'C.ht']:
    if os.path.exists(file):
        shutil.rmtree(file)
