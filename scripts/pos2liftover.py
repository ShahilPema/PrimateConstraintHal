import hail as hl
import argparse
import shutil
import os
import glob

parser = argparse.ArgumentParser(description="Perform liftover from species-specific BED to human BED using Hail.")
parser.add_argument("--species", required=True, help="Name of the species to process (e.g., Pithecia_pithecia).")
parser.add_argument("--species_fasta_path", required=True, help="Path to the species-specific reference FASTA file.")
parser.add_argument("--species_index_path", required=True, help="Path to the index file for the species-specific reference FASTA.")
parser.add_argument("--pos_ht", required=True, help="Path to the regionsxvar Hail table (regionsxvar.ht).")
parser.add_argument("--vep_annotations", required=True, help="Path to the vep annotations Hail table (vep.ht).")
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

complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}

lifttable = hl.read_table(args.pos_ht)

lifttable = lifttable.select(
    locus = lifttable.locus,
    species_locus = lifttable[f'locus_{args.species}'],
    hg38_chr = lifttable.locus.contig, 
    hg38_position = lifttable.locus.position,
    species_contig = lifttable[f'locus_{args.species}'].contig,
    species_position = lifttable[f'locus_{args.species}'].position,
    species_pentamer = lifttable[f'{args.species}_pos_info'].species_pentamer,
    strand_v_human_ref = lifttable[f'{args.species}_pos_info'].strand_v_human_ref,
    hg38_ref = lifttable[f'{args.species}_pos_info'].hg38_pentamer[2],
    species_ref = lifttable[f'{args.species}_pos_info'].species_pentamer[2],
    ref_match = lifttable[f'{args.species}_pos_info'].ref_match,
    overlapping_annotations = lifttable[f'{args.species}_pos_info'].overlapping_annotations
)

lifttable = lifttable.checkpoint('cpt3_0.ht', overwrite=True)

lifttable = lifttable.key_by('hg38_chr', 'hg38_position')

lifttable = lifttable.checkpoint('cpt3_1.ht', overwrite=True)

vep = hl.read_table(args.vep_annotations)
vep = vep.key_by()
vep = vep.annotate(hg38_position = hl.int32(vep.hg38_position))
vep = vep.key_by('hg38_chr', 'hg38_position')

vep = vep.checkpoint('cpt3_2.ht', overwrite=True)

lifttable = lifttable.join(vep, how='left')

lifttable = lifttable.add_index()

lifttable = lifttable.checkpoint('cpt3_3.ht', overwrite=True)

codon_groups = lifttable.filter(lifttable.region_type == 'CDS')

codon_groups = codon_groups.checkpoint('cpt3_4.ht', overwrite=True)

# Group by transcript, region_type, and protein_start
codon_groups = codon_groups.group_by(
    'transcript', 'region_type', 'protein_start'
).aggregate(
    species_loci=hl.agg.collect(hl.struct(
        contig=codon_groups.species_contig,
        position=codon_groups.species_position,
        codon_position=codon_groups.codon_position,
        strand_v_human_transcript=codon_groups.strand_v_human_ref * codon_groups.strand,
        idx=codon_groups.idx
    ))
)

codon_groups = codon_groups.checkpoint('cpt4.ht', overwrite=True)

# Sort loci by contig and position
codon_groups = codon_groups.annotate(
    loci_by_contig_sorted=hl.array([hl.sorted(codon_groups.species_loci, key=lambda x: x.position * x.strand_v_human_transcript)])
)

# Process codon windows
codon_groups = codon_groups.annotate(
    processed_rows=codon_groups.loci_by_contig_sorted.flatmap(lambda sorted_loci:
        hl.range(hl.len(sorted_loci) - 2).map(lambda i:
            hl.struct(
                contig=sorted_loci[i].contig,
                position=sorted_loci[i].position,
                idx=sorted_loci[i].idx,
                window=hl.if_else(
                    (
                        #A valid codon has 3 positions in the correct order on the same strand
                        (i <= hl.len(sorted_loci) - 3) & 
                        (hl.set([[0, 1, 2], [2, 1, 0]]).contains(hl.range(3).map(lambda j: sorted_loci[i+j].codon_position))) &
                        (hl.len(hl.set(hl.range(3).map(lambda j: sorted_loci[i+j].strand_v_human_transcript))) == 1)
                    ), hl.if_else(
                        #A valid codon must have at least 2 adjacent positions
                        (
                            (hl.abs(sorted_loci[i].position - sorted_loci[i+1].position) == 1) |
                            (hl.abs(sorted_loci[i+1].position - sorted_loci[i+2].position) == 1)
                        ),
                        sorted_loci[i:i+3],
                        hl.missing(hl.tarray(sorted_loci[i].dtype))
                    ),
                    #codons that fail will result in a missing array and be annotated as 'incomplete' or 'nonhomo' if position is missing
                    hl.missing(hl.tarray(sorted_loci[i].dtype))
                )
            )
        )
    )
)

codon_groups = codon_groups.checkpoint('cpt5.ht', overwrite=True)

struct_fields = dict(codon_groups.processed_rows.first().dtype.items())
w_added_fields = {**struct_fields, **{'idxs': hl.tarray(hl.tint32), 'species_codon_flag': hl.tstr, 'species_codon': hl.tstr}}  # Unpack both dicts
new_missing_placeholder = hl.tstruct(**w_added_fields)

species_rg = hl.ReferenceGenome.from_fasta_file(
    name=args.species,
    fasta_file=args.species_fasta_path,
    index_file=args.species_index_path
)

# Annotate codon flags and sequences
codon_groups = codon_groups.annotate(
    processed_rows=codon_groups.processed_rows.map(lambda struct:
        hl.if_else(
            hl.is_defined(struct.window),
            struct.annotate(
                idxs = hl.array([struct.window[0].idx, struct.window[1].idx, struct.window[2].idx]),
                species_codon_flag=hl.case()
                    .when((
                            (hl.abs(struct.window[1].position - struct.window[0].position) > 10000)
                            | (hl.abs(struct.window[2].position - struct.window[1].position) > 10000)
                    ), "discontiguous_largegap")
                    .when((
                            #previous code already ensures that at least 2 adjacent positions exist
                            (hl.abs(struct.window[1].position - struct.window[0].position) != 1)
                            | (hl.abs(struct.window[2].position - struct.window[1].position) != 1)
                    ), "discontiguous")
                    .default("contiguous"),
                species_codon=hl.if_else(
                    struct.window.first().strand_v_human_transcript == -1,
                    hl.delimit(struct.window.map(lambda x: hl.get_sequence(x.contig, x.position, reference_genome=args.species).upper()), '').translate(complement),
                    hl.delimit(struct.window.map(lambda x: hl.get_sequence(x.contig, x.position, reference_genome=args.species).upper()), '')
                )
            ),
            hl.missing(new_missing_placeholder)
        )
    )
)

codon_groups = codon_groups.key_by()

codon_groups = codon_groups.select('processed_rows')

codon_groups = codon_groups.explode('processed_rows')

codon_groups = codon_groups.annotate(
    idxs = codon_groups.processed_rows.idxs,
    species_codon = codon_groups.processed_rows.species_codon,
    species_codon_flag = codon_groups.processed_rows.species_codon_flag
).drop('processed_rows')

codon_groups = codon_groups.checkpoint('cpt6.ht', overwrite=True)

# Create a table for valid codons instead of a dictionary
codon_groups = codon_groups.explode('idxs')

codon_groups = codon_groups.rename({'idxs': 'idx'})

codon_groups = codon_groups.filter(
    hl.set(["contiguous", "discontiguous", "discontiguous_largegap"]).contains(codon_groups.species_codon_flag)
)

codon_groups = codon_groups.key_by('idx')

codon_groups = codon_groups.checkpoint('cpt7.ht', overwrite=True)

# Final annotated table
lifttable = lifttable.annotate(
    species_codon_flag=codon_groups[lifttable.idx].species_codon_flag,
    species_codon=codon_groups[lifttable.idx].species_codon
)

lifttable = lifttable.annotate(
    species_codon_flag = hl.if_else(
        hl.is_missing(lifttable.species_position) & (lifttable.region_type == 'CDS'),
        'nonhomo',
        lifttable.species_codon_flag
    ),
    species_codon = hl.if_else(
        hl.is_missing(lifttable.species_position) & (lifttable.region_type == 'CDS'),
        hl.missing(hl.tstr),
        lifttable.species_codon
    )
)

lifttable = lifttable.annotate(
    species_codon_flag = hl.if_else(
        hl.is_missing(lifttable.species_codon_flag) & (lifttable.region_type == 'CDS'),
        'incomplete',
        lifttable.species_codon_flag
    ),
    species_codon = hl.if_else(
        hl.is_missing(lifttable.species_codon) & (lifttable.region_type == 'CDS'),
        hl.missing(hl.tstr),
        lifttable.species_codon
    ),
    human_codon = lifttable.codons.upper()
)

lifttable = lifttable.key_by()

lifttable = lifttable.drop('idx', 'hg38_chr', 'hg38_position', 'codons')

lifttable = lifttable.checkpoint('cpt8.ht', overwrite=True)

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

lifttable = lifttable.checkpoint('cpt9.ht', overwrite=True)

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

def annotate_alt_consequence(ref_aa, alt_aa, codon_position):
    return hl.case() \
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
    species_alt_cons=annotate_alt_consequence(lifttable.species_aa, lifttable.species_alt_aa, lifttable.codon_position),
    codon_mismatch_cons = hl.if_else(
        lifttable.codon_match == 0,
        annotate_alt_consequence(lifttable.human_aa, lifttable.species_aa, lifttable.codon_position),
        hl.missing(hl.tstr)
    )
)

lifttable = lifttable.annotate(
    alt_cons_match = hl.case() \
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
              .when(lifttable.strand_v_human_ref == 1,lifttable.alleles) \
              .when(lifttable.strand_v_human_ref == -1, hl.map(lambda x: x.translate(complement), lifttable.alleles)) \
              .default(hl.missing(hl.tarray(hl.tstr)))
             )
        .when(lifttable.ref_match == 0,
              hl.case() \
              .when(lifttable.strand_v_human_ref == 1, process_alleles(lifttable.alleles, lifttable.species_ref)) \
              .when(lifttable.strand_v_human_ref == -1, process_alleles(hl.map(lambda x: x.translate(complement), lifttable.alleles), lifttable.species_ref)) \
              .default(hl.missing(hl.tarray(hl.tstr)))
             )
        .default(hl.missing(hl.tarray(hl.tstr))),
)

lifttable = lifttable.checkpoint('cpt10.ht', overwrite=True)

lifttable_no_overlap = lifttable.filter(lifttable.overlapping_annotations == 1)

lifttable_no_overlap = lifttable_no_overlap.annotate(
    transcript_mapping = hl.set([hl.struct(
            transcript = lifttable_no_overlap.transcript,
            region_type = lifttable_no_overlap.region_type,
            strand = lifttable_no_overlap.strand,
            human_codon = lifttable_no_overlap.human_codon,
            species_codon = lifttable_no_overlap.species_codon,
            species_codon_flag = lifttable_no_overlap.species_codon_flag,
            codon_match = lifttable_no_overlap.codon_match,
            aa_match = lifttable_no_overlap.aa_match,
            codon_mismatch_cons = lifttable_no_overlap.codon_mismatch_cons,
            human_alt_cons = lifttable_no_overlap.human_alt_cons,
            species_alt_cons = lifttable_no_overlap.species_alt_cons,
            alt_cons_match = lifttable_no_overlap.alt_cons_match
        )])
)

lifttable_no_overlap = lifttable_no_overlap.select('locus', 'alleles', 'species_locus', 'species_alleles', 'transcript_mapping')


lifttable_overlap = lifttable.filter(lifttable.overlapping_annotations > 1)
lifttable_overlap = lifttable_overlap.group_by('locus', 'alleles', 'species_locus', 'species_alleles').aggregate(
    transcript_mapping=hl.agg.collect_as_set(
        hl.struct(
            transcript = lifttable_overlap.transcript,
            region_type = lifttable_overlap.region_type,
            strand = lifttable_overlap.strand,
            human_codon = lifttable_overlap.human_codon,
            species_codon = lifttable_overlap.species_codon,
            species_codon_flag = lifttable_overlap.species_codon_flag,
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

lifttable = lifttable.rename({ 'species_locus' :f'{args.species}_locus', 'species_alleles' : f'{args.species}_alleles', 'transcript_mapping' : f'{args.species}_transcript_mapping'})

lifttable = lifttable.repartition(int(args.cpus))

lifttable = lifttable.key_by(f'{args.species}_locus', f'{args.species}_alleles')

lifttable.write(f'{args.output}/{args.species}_var_info.ht', overwrite=True)


files_to_remove = glob.glob('cpt*.ht')
for file in files_to_remove:
    if os.path.exists(file):
        shutil.rmtree(file)