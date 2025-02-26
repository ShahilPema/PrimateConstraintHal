import hail as hl
import psutil
import argparse
import shutil
import os
import math

parser = argparse.ArgumentParser(description="Aggregate counts across all references.")
parser.add_argument("--vep", required=True, help="VEP hail table, containing all sites")
parser.add_argument("--human_fa", required=True, help="Human reference fasta file")
parser.add_argument("--human_faidx", required=True, help="Human reference fasta index file")
parser.add_argument("--var_hts", required=True, help="List of aggregated varaiant-level matrix tables for each reference species")
parser.add_argument("--pos_hts", required=True, help="List of aggregated position-level matrix tables for each reference species")
parser.add_argument("--mapping", required=True, help="233primates_s1.csv")
args = parser.parse_args()

#NOTE change cpus so each partition is more than 128MB. 

hl_cpus = 32

forks = 1

memory = int(psutil.virtual_memory().total/(1024 ** 3)*0.85/forks)

config = {
    'spark.driver.memory': f'{memory}g',  #Set to total memory
    'spark.executor.memory': f'{memory}g'  #Set to total memory
}

hl.init(spark_conf=config, master=f'local[32]')

id_info = hl.import_table(args.mapping, delimiter=',')
id_info = id_info.rename({'\ufeffID': 'ID'})
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
order_dict = hl.dict(hl.zip(id_info.ID.collect(),
                          [family_to_order.get(x) for x in id_info.FAMILY.collect()]))

human_rg = hl.ReferenceGenome.from_fasta_file(
    name='Homo_sapiens',
    fasta_file='Homo_sapiens.fa',
    index_file='Homo_sapiens.fa.fai'
)



all_positions = hl.read_table(args.vep)
all_positions = all_positions.filter(all_positions.locus.contig == 'chr21')
all_positions = all_positions.select('regions')
all_positions = all_positions.key_by()
all_positions = all_positions.annotate(
    locus = hl.locus(all_positions.locus.contig, all_positions.locus.position, reference_genome='Homo_sapiens')
)
all_positions = all_positions.key_by('locus', 'alleles')


var_level_files = args.var_hts.split()

species_count_dict = {}
species_list = []

for file in var_level_files:
    ht = hl.read_table(file)

    ht = ht.distinct()

    tcall_columns = [name for name, dtype in ht.row.dtype.items() if dtype == hl.tcall]

    species = file.split('/')[-1].removesuffix('_lifted_vars.ht')
    species_list.append(species)
    species_count_dict[species] = len(tcall_columns)
    
    #Apply hl.agg.call_stats to this combined array
    ht = ht.annotate(
        species_data=ht[f'{species}_data'].annotate(
            counts=hl.if_else(
                hl.is_defined(ht.alleles),
                hl.struct(
                    AC=hl.sum([ht[col].n_alt_alleles() for col in tcall_columns]),
                    AN=hl.sum([hl.is_defined(ht[col]) for col in tcall_columns]) * 2,
                    HZ=hl.sum([ht[col].is_hom_var() for col in tcall_columns]),
                    HET=hl.sum([ht[col].is_het() for col in tcall_columns]),
                    samples = hl.array(
                                hl.filter(
                                    lambda x: hl.is_defined(x),  # Keep only defined (non-missing) values
                                    [hl.if_else(ht[col].is_non_ref(), col, hl.missing(hl.tstr)) for col in tcall_columns]
                                )
                            ),
                    # Add sample IDs for homozygotes and heterozygotes
                    het_samples = hl.array(
                                hl.filter(
                                    lambda x: hl.is_defined(x),
                                    [hl.if_else(ht[col].is_het(), col, hl.missing(hl.tstr)) for col in tcall_columns]
                                )
                            ),
                    hom_samples = hl.array(
                                hl.filter(
                                    lambda x: hl.is_defined(x),
                                    [hl.if_else(ht[col].is_hom_var(), col, hl.missing(hl.tstr)) for col in tcall_columns]
                                )
                            ),
                ),
                hl.missing(hl.tstruct(
                    AC=hl.tint, 
                    AN=hl.tint, 
                    HZ=hl.tint, 
                    HET=hl.tint,  # Add heterozygote type
                    samples=hl.tarray(hl.tstr),
                    het_samples=hl.tarray(hl.tstr),
                    hom_samples=hl.tarray(hl.tstr)
                ))
            )
        )
    )

    ht = ht.drop(*tcall_columns, f'{species}_data')

    ht = ht.rename({'species_data': f'{species}_data'})

    ht = ht.checkpoint(f'./tmp/{species}_data.ht', overwrite=True)
    

for species in species_list:
    ht = hl.read_table(f'./tmp/{species}_data.ht')
    all_positions = all_positions.annotate(
        species = ht[all_positions.key][f'{species}_data']
    )
    all_positions = all_positions.rename({'species': f'{species}_data'})
    
    all_positions = all_positions.checkpoint(f'./tmp/{species}_data_2.ht', overwrite=True)

all_positions = all_positions.key_by('locus')

all_positions = all_positions.checkpoint('./tmp/checkpoint1.ht', overwrite=True)

for species in species_list:

    ht = hl.read_table(f'./tmp/{species}_data.ht')

    ht = ht.filter(hl.is_defined(ht[f'{species}_data'].species_locus) & hl.is_missing(ht[f'{species}_data'].species_alleles))

    ht = ht.key_by('locus')

    all_positions = all_positions.annotate(
        species = hl.if_else(
            hl.is_missing(all_positions[f'{species}_data'].species_locus),
            ht[all_positions.key][f'{species}_data'],
            all_positions[f'{species}_data']
        )
    )

    all_positions = all_positions.annotate(
        species = all_positions.species.annotate(
            transcript_mappings = hl.if_else(
                hl.is_missing(all_positions[f'{species}_data'].species_alleles),
                hl.missing('set<struct{transcript: str, region_type: str, strand: int32, codon_match: int32, aa_match: int32, codon_mismatch_cons: str, human_alt_cons: str, species_alt_cons: str, alt_cons_match: int32}>'),
                all_positions[f'{species}_data'].transcript_mappings
            )
        )
    )
    
    all_positions = all_positions.annotate(
        species = all_positions.species.annotate(
            transcript_mappings = all_positions.species.transcript_mappings.filter(
                lambda x: x.region_type == 'CDS'
            )
        )
    )
    
    all_positions = all_positions.annotate(
        species = all_positions.species.annotate(
            transcript_mappings = hl.if_else(
                hl.is_missing(all_positions.species.transcript_mappings),
                hl.missing('set<struct{transcript: str, region_type: str, strand: int32, codon_match: int32, aa_match: int32, codon_mismatch_cons: str, human_alt_cons: str, species_alt_cons: str, alt_cons_match: int32, annotation: str}>'),
                all_positions.Aotus_nancymaae_data.transcript_mappings.map(
                    lambda x: x.annotate(annotation=x.transcript + '_' + x.region_type)
                )
            )
        )
    )
    
    
            
    all_positions = all_positions.drop(f'{species}_data')

    all_positions = all_positions.rename({'species': f'{species}_data'})   
    
    all_positions = all_positions.checkpoint(f'./tmp/{species}_data_2.ht', overwrite=True)


all_positions = hl.read_table('./tmp/Sapajus_apella_data_2.ht')
species_list = list(species_count_dict.keys())

species_count_dict = hl.literal(species_count_dict)



# Compute species-level summaries at the row level
all_positions = all_positions.annotate(
    species_locus_defined=hl.sum([
        hl.is_defined(all_positions[f'{species}_data'].species_locus) for species in species_list
    ]),
    species_alleles_defined=hl.sum([
        hl.is_defined(all_positions[f'{species}_data'].species_alleles) for species in species_list
    ])
)

all_positions = all_positions.annotate(
    sumAC=hl.sum([
        hl.or_else(all_positions[f'{species}_data'].counts.AC, 0) for species in species_list
    ]),
    sumAN=hl.sum([
        hl.or_else(all_positions[f'{species}_data'].counts.AN, 0) for species in species_list
    ]),
    sumHZ=hl.sum([
        hl.or_else(all_positions[f'{species}_data'].counts.HZ, 0) for species in species_list
    ]),
    sumHET=hl.sum([  # Add total heterozygote count
        hl.or_else(all_positions[f'{species}_data'].counts.HET, 0) for species in species_list
    ]),
    samples = hl.flatten([
        hl.or_else(all_positions[f'{species}_data'].counts.samples, hl.empty_array(hl.tstr)) for species in species_list
    ]),
    # Add flattened lists of heterozygous and homozygous samples
    het_samples = hl.flatten([
        hl.or_else(all_positions[f'{species}_data'].counts.het_samples, hl.empty_array(hl.tstr)) for species in species_list
    ]),
    hom_samples = hl.flatten([
        hl.or_else(all_positions[f'{species}_data'].counts.hom_samples, hl.empty_array(hl.tstr)) for species in species_list
    ]),
    sumAC_ref_match=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 1),
            hl.or_else(all_positions[f'{species}_data'].counts.AC, 0),
            0
        ) for species in species_list
    ]),
    sumAN_ref_match=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 1),
            hl.or_else(all_positions[f'{species}_data'].counts.AN, 0),
            0
        ) for species in species_list
    ]),
    sumHZ_ref_match=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 1),
            hl.or_else(all_positions[f'{species}_data'].counts.HZ, 0),
            0
        ) for species in species_list
    ]),
    sumHET_ref_match=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 1),
            hl.or_else(all_positions[f'{species}_data'].counts.HET, 0),
            0
        ) for species in species_list
    ]),
    sumAC_ref_mismatch=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 0),
            hl.or_else(all_positions[f'{species}_data'].counts.AC, 0),
            0
        ) for species in species_list
    ]),
    sumAN_ref_mismatch=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 0),
            hl.or_else(all_positions[f'{species}_data'].counts.AN, 0),
            0
        ) for species in species_list
    ]),
    sumHZ_ref_mismatch=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 0),
            hl.or_else(all_positions[f'{species}_data'].counts.HZ, 0),
            0
        ) for species in species_list
    ]),
    sumHET_ref_mismatch=hl.sum([
        hl.if_else(
            hl.is_defined(all_positions[f'{species}_pos_info']) & (all_positions[f'{species}_pos_info'].ref_match == 0),
            hl.or_else(all_positions[f'{species}_data'].counts.HET, 0),
            0
        ) for species in species_list
    ])
)

all_positions = all_positions.checkpoint('all_positions_cpt2.ht', overwrite=True)

all_positions = all_positions.annotate(
    AF=(all_positions.sumAC / all_positions.sumAN),
    species=hl.array(hl.set(all_positions.samples.map(lambda x: species_dict[x]))),
    family=hl.array(hl.set(all_positions.samples.map(lambda x: family_dict[x]))),
    reference=hl.array(hl.set(all_positions.samples.map(lambda x: reference_dict[x]))),
    # Add order annotation
    order=hl.array(hl.set(all_positions.samples.map(lambda x: order_dict[x])))
)

all_positions = all_positions.checkpoint('all_positions_cpt3.ht', overwrite=True)

# Create order-specific homozygote and heterozygote counts
orders = ["Platyrrhini", "Catarrhini", "Strepsirrhini", "Tarsiiformes"]

# Function to count samples by order
def count_by_order(sample_list, order_name):
    return hl.sum(sample_list.map(lambda x: hl.if_else(order_dict[x] == order_name, 1, 0)))

# Add order-specific counts
all_positions = all_positions.annotate(
    order_counts=hl.struct(**{
        f"{order}_hom": count_by_order(all_positions.hom_samples, order) for order in orders
    }, **{
        f"{order}_het": count_by_order(all_positions.het_samples, order) for order in orders
    })
)

# Add var_primate_orders column - array of orders with variants
all_positions = all_positions.annotate(
    var_primate_orders=hl.array(hl.set(
        all_positions.het_samples.map(lambda x: order_dict[x]).extend(
            all_positions.hom_samples.map(lambda x: order_dict[x])
        )
    ))
)


fields = ['codon_match', 'aa_match', 'alt_cons_match']
all_fields = ['annotation'] + fields

# Combine all arrays, dynamically using collected species names
all_positions = all_positions.annotate(
    combined=hl.flatten([
        hl.or_else(
            hl.array(all_positions[f"{species}_data"].transcript_mappings.map(lambda x: x.select(*all_fields).annotate(original_column=species))),
            hl.empty_array(hl.tstruct(annotation=hl.tstr, **{field: hl.tint32 for field in fields}, original_column=hl.tstr))
        )
        for species in species_list
    ])
)

# Annotate the table with groups dynamically for each row
all_positions = all_positions.annotate(
    groups=hl.set(all_positions.combined.map(lambda x: x.annotation))
)

all_positions = all_positions.checkpoint('./tmp/all_positions_cpt4.ht', overwrite=True)

# Calculate the sum and weighted mean for each group
def calculate_struct(group):
    used_columns = hl.set(all_positions.combined.filter(lambda x: x.annotation == group).map(lambda x: x.original_column))
    denominator = hl.sum(used_columns.map(lambda col: species_count_dict[col]))
    return hl.struct(
        annotation=group,
        **{
            f"{f}_def_species": hl.sum(all_positions.combined.filter(lambda x: x.annotation == group).map(lambda x: x[f])) for f in fields
        },
        **{
            f"{f}_weighted_mean": (
                hl.sum(all_positions.combined.filter(lambda x: x.annotation == group).map(lambda x: x[f] * species_count_dict[x.original_column])) / denominator
            ) for f in fields
        }
    )

all_positions = all_positions.annotate(
    transcript_mappings=hl.array(all_positions.groups.map(calculate_struct))
)

all_positions = all_positions.drop('combined', 'groups', *[f"{species}_data" for species in species_list])

all_positions.checkpoint('./tmp/aggregated_primate_counts.ht', overwrite=True)

pos_level_files = args.pos_hts.split()

# Step 2: Annotate all_positions with position-level data
for file in pos_level_files:
    # Extract species name from the file name
    species = file.split('/')[-1].removesuffix('_pos_info.ht')
    
    # Read the position-level data
    pos_ht = hl.read_table(file)
    
    # Ensure the data is keyed by Homo_sapiens locus for proper join
    pos_ht = pos_ht.key_by('locus')
    
    # Annotate the all_positions table with the position-level data
    all_positions = all_positions.annotate(
        **{f'{species}_pos_info': pos_ht[all_positions.key][f'{species}_pos_info']}
    )
    
    all_positions = all_positions.checkpoint(f'./tmp/{species}_data.ht', overwrite=True)

# Step 3: Calculate aggregated statistics
fields = ['ref_match', 'min_hamming', 'pentamer_error']
all_positions = all_positions.annotate(
    **{
        f'sum_{field}': hl.sum([
            hl.or_else(all_positions[f'{species}_pos_info'][field], 0) for species in species_list
        ])
        for field in ['ref_match']
    },
    **{
        f'mean_{field}': (
            hl.sum([
                hl.or_else(all_positions[f'{species}_pos_info'][field], 0) for species in species_list
            ]) / 
            hl.sum([
                hl.if_else(hl.is_defined(all_positions[f'{species}_pos_info'][field]), 1, 0) for species in species_list
            ])
        )
        for field in ['min_hamming', 'pentamer_error']
    },
    **{
        f'weighted_mean_{field}': (
            hl.sum([
                hl.or_else(all_positions[f'{species}_pos_info'][field], 0) * species_count_dict[species]
                for species in species_list
            ]) / 
            hl.sum([
                hl.if_else(hl.is_defined(all_positions[f'{species}_pos_info'][field]), species_count_dict[species], 0) for species in species_list
            ])
        )
        for field in ['min_hamming', 'pentamer_error']
    }
)

all_positions = all_positions.drop(*[f'{species}_pos_info' for species in species_list] + ['samples'])

# Step 5: Save the updated table
all_positions.write('aggregated_primate_counts_with_position_info.ht', overwrite=True)

shutil.rmtree('tmp')
