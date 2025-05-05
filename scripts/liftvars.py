import hail as hl
import argparse
import shutil
import os

parser = argparse.ArgumentParser(description="Perform liftover from species-specific BED to human BED using Hail.")
parser.add_argument("--species", required=True, help="Name of the species to process (e.g., Pithecia_pithecia).")
parser.add_argument("--sample_mts", required=True, help="Sample hail tables to be merged ")
parser.add_argument("--species_ht", required=True, help="Species variant info hail table.")
parser.add_argument("--pos_ht", required=True, help="Position-level variant info hail table.")
parser.add_argument("--mapping", required=True, help="233primates_s1.csv")
parser.add_argument("--cpus", required=True, help="cpus for task")
parser.add_argument("--memory", required=True, help="memory for task")
parser.add_argument("--contig_count", required=True, help="Total contigs.")
parser.add_argument("--tmpdir", required=True, help="Temporary directory")
parser.add_argument("--outdir", required=True, help="Total system cpus")
args = parser.parse_args()

config = {
    'spark.driver.memory': f'{args.memory}',  #Set to total memory
    'spark.executor.memory': f'{args.memory}',
    'spark.local.dir': args.tmpdir,
    'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}',
    'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}'
}

hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)

id_info = hl.import_table(args.mapping, delimiter=',', impute=True)
id_info = id_info.rename({'???ID': 'ID'})
species_dict = hl.dict(hl.zip(id_info.ID.collect(), id_info.SPECIES.collect()))
family_dict = hl.dict(hl.zip(id_info.ID.collect(), id_info.FAMILY.collect()))
reference_dict = hl.dict(hl.zip(id_info.ID.collect(), id_info.REFERENCE.collect()))

# Add family to order mapping
family_to_order = hl.literal({
    "Pitheciidae": "Platyrrhini",  # New World Monkeys
    "Callitrichidae": "Platyrrhini",
    "Atelidae": "Platyrrhini",
    "Cebidae": "Platyrrhini",
    "Aotidae": "Platyrrhini",
    "Cercopithecidae": "Catarrhini",  # Old World Monkeys
    "Hylobatidae": "Catarrhini",  # Lesser Apes
    "Hominidae": "Catarrhini",  # Great Apes
    "Lorisidae": "Strepsirrhini",  # Lorises
    "Indriidae": "Strepsirrhini",  # Lemurs
    "Daubentoniidae": "Strepsirrhini",
    "Lemuridae": "Strepsirrhini",
    "Galagidae": "Strepsirrhini",  # Galagos
    "Lepilemuridae": "Strepsirrhini",
    "Cheirogaleidae": "Strepsirrhini",
    "Tarsiidae": "Tarsiiformes"  # Tarsiers
})

# Create order dictionary based on family
order_dict = hl.dict(hl.zip(
                            id_info.ID.collect(),
                            [family_to_order.get(x) for x in id_info.FAMILY.collect()]
                        )
                    )

species_ht = hl.read_table(args.species_ht)

species_ht = species_ht.annotate(
    contig = species_ht[f'{args.species}_locus'].contig,
    position = species_ht[f'{args.species}_locus'].position
)

contig_count = int(args.contig_count)
if contig_count > 100000: 
	species_ht = species_ht.key_by('contig', 'position', f'{args.species}_alleles')

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
    species_ht = species_ht.annotate(GT = ht[species_ht.key].GT)
    species_ht = species_ht.annotate(GT = hl.case() \
                                .when(hl.is_missing(species_ht.alleles), hl.missing(hl.tcall)) \
                                .when(hl.is_defined(species_ht.GT), species_ht.GT) \
                                .default(hl.call(0,0))
    )
    species_ht = species_ht.rename({'GT':sample_id})
    species_ht = species_ht.checkpoint(f'{args.tmpdir}/checkpoint_{count}.ht', overwrite=True)
    count += 1

tcall_columns = [name for name, dtype in species_ht.row.dtype.items() if dtype == hl.tcall]

species_ht = species_ht.annotate(
    species_locus = species_ht[f'{args.species}_locus'],
    species_alleles = species_ht[f'{args.species}_alleles'],
    transcript_mappings = species_ht[f'{args.species}_transcript_mapping'],
    samples = hl.if_else(
        hl.is_defined(species_ht.alleles),
        hl.array(
            hl.filter(
                lambda x: hl.is_defined(x),  # Keep only defined (non-missing) values
                [hl.if_else(species_ht[col].is_non_ref(), col, hl.missing(hl.tstr)) for col in tcall_columns]
            )
        ),
        hl.missing(hl.tarray(hl.tstr))
    ),
    het_samples = hl.if_else(
        hl.is_defined(species_ht.alleles),
        hl.array(
            hl.filter(
                lambda x: hl.is_defined(x),  # Keep only defined (non-missing) values
                [hl.if_else(species_ht[col].is_het(), col, hl.missing(hl.tstr)) for col in tcall_columns]
            )
        ),
        hl.missing(hl.tarray(hl.tstr))
    ),
    hom_samples = hl.if_else(
        hl.is_defined(species_ht.alleles),
        hl.array(
            hl.filter(
                lambda x: hl.is_defined(x),
                [hl.if_else(species_ht[col].is_hom_var(), col, hl.missing(hl.tstr)) for col in tcall_columns]
            )
        ),
        hl.missing(hl.tarray(hl.tstr))
    ),
    ref_samples = hl.if_else(
        hl.is_defined(species_ht.alleles),
        hl.array(
            hl.filter(
                lambda x: hl.is_defined(x), 
                [hl.if_else(species_ht[col].is_hom_ref(), col, hl.missing(hl.tstr)) for col in tcall_columns]
            )
        ),
        hl.missing(hl.tarray(hl.tstr))
    )
)

species_ht = species_ht.annotate(
    species = hl.set(species_ht.samples.map(lambda x: species_dict[x])),
    family = hl.set(species_ht.samples.map(lambda x: family_dict[x])),
    reference = hl.set(species_ht.samples.map(lambda x: reference_dict[x])),
    order = hl.set(species_ht.samples.map(lambda x: order_dict[x])),
    n_hom = hl.or_else(hl.len(species_ht.hom_samples), 0),
    n_het = hl.or_else(hl.len(species_ht.het_samples), 0),
    n_ref = hl.or_else(hl.len(species_ht.ref_samples), 0)
    )

species_ht = species_ht.annotate(
    AC = species_ht.n_hom * 2 + species_ht.n_het,
    AN = species_ht.n_hom * 2 + species_ht.n_het * 2 + species_ht.n_ref * 2,
)


species_ht = species_ht.checkpoint(f'{args.tmpdir}/{args.species}_data1.ht', overwrite=True)

species_ht = species_ht.annotate(
    AF = species_ht.AC / species_ht.AN
)

orders = ["Platyrrhini", "Catarrhini", "Strepsirrhini", "Tarsiiformes"]

def count_by_order(sample_list, order_name):
    return hl.sum(sample_list.map(lambda x: hl.if_else(order_dict[x] == order_name, 1, 0)))

# Add order-specific counts
species_ht = species_ht.annotate(
    **{
        f"{order}_hom": count_by_order(hl.array(hl.set(species_ht.hom_samples)), order) for order in orders
    }, **{
        f"{order}_het": count_by_order(hl.array(hl.set(species_ht.het_samples)), order) for order in orders
    }
)

species_ht = species_ht.checkpoint(f'{args.tmpdir}/{args.species}_data1_2.ht', overwrite=True)

species_ht = species_ht.key_by('species_locus')

# Read the position-level data
pos_ht = hl.read_table(args.pos_ht)
    
# Ensure the data is keyed by Homo_sapiens locus for proper join
pos_ht = pos_ht.key_by(f'locus_{args.species}')

species_ht = species_ht.checkpoint(f'{args.tmpdir}/{args.species}_data2.ht', overwrite=True)

 
# Annotate the all_positions table with the position-level data
species_ht = species_ht.annotate(
    pos_info = pos_ht[species_ht.key][f'{args.species}_pos_info']
)

species_ht = species_ht.checkpoint(f'{args.tmpdir}/{args.species}_data2_1.ht', overwrite=True)
        
# Annotate with position-level info 'ref_match', 'min_hamming', 'pentamer_error' and 'pentamer'
species_ht = species_ht.annotate(
    ref_match = species_ht.pos_info.ref_match,
    min_hamming = species_ht.pos_info.min_hamming,
    pentamer_error = species_ht.pos_info.pentamer_error,
    species_pentamer = species_ht.pos_info.species_pentamer,
    reciprocal_mapping = species_ht.pos_info.reciprocal_mapping,
    human2species_mappings = species_ht.pos_info.human2species_mappings,
    overlapping_annotations = species_ht.pos_info.overlapping_annotations,
    strand_v_human_ref = species_ht.pos_info.strand_v_human_ref
)

species_ht = species_ht.checkpoint(f'{args.tmpdir}/{args.species}_data.ht', overwrite=True)
    
non_call_columns = [
    'species_locus', 'species_alleles', 'transcript_mappings', 'species', 'family', 'order', 'reference', 'AN', 'AC', 
    'AF', 'n_hom', 'n_het', 'samples', 'hom_samples', 'het_samples', 'ref_match', 'min_hamming', 'pentamer_error', 
    'species_pentamer', 'reciprocal_mapping', 'human2species_mappings', 'overlapping_annotations', 'strand_v_human_ref'
    ] + [f"{order}_hom" for order in orders] + [f"{order}_het" for order in orders]

species_ht = species_ht.annotate(
    species_data = hl.struct(
        **{col: species_ht[col] for col in non_call_columns}
    ),
    call_data = hl.struct(
        **{col: species_ht[col] for col in tcall_columns}
    )
)
species_ht = species_ht.key_by()

species_ht = species_ht.select('locus', 'alleles', 'species_data', 'call_data')

species_ht = species_ht.checkpoint(f'{args.tmpdir}/A.ht', overwrite=True)

species_ht = species_ht.group_by('locus', 'alleles').aggregate(
    species_data = hl.agg.collect(species_ht.species_data),
    call_data = hl.agg.collect(species_ht.call_data)
)

species_ht = species_ht.rename({'species_data': f'{args.species}_data'})

species_ht.write(f'{args.outdir}/{args.species}_lifted_vars.ht', overwrite=True)

column_data = column_data.repartition(1)

column_data.write(f'{args.outdir}/{args.species}_colinfo_vars.ht', overwrite=True)
