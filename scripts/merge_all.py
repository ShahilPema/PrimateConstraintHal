import hail as hl
import argparse

parser = argparse.ArgumentParser(description="Aggregate counts across all references.")
parser.add_argument("--vep", required=True, help="VEP hail table, containing all sites")
parser.add_argument("--human_fa", required=True, help="Human reference fasta file")
parser.add_argument("--human_faidx", required=True, help="Human reference fasta index file")
parser.add_argument("--var_hts", required=True, help="List of aggregated varaiant-level matrix tables for each reference species")
parser.add_argument("--cpus", required=True, help="Number of CPUs to use")
parser.add_argument("--memory", required=True, help="Memory to use")
parser.add_argument("--tmpdir", required=True, help="Temporary directory")
parser.add_argument("--outdir", required=True, help="Output directory")
args = parser.parse_args()

config = {
    'spark.driver.memory': f'{args.memory}',  #Set to total memory
    'spark.executor.memory': f'{args.memory}',  #Set to total memory
    'spark.local.dir': args.tmpdir
}

hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)


all_positions = hl.read_table(args.vep)
all_positions = all_positions.select('regions')

var_level_files = args.var_hts.split()

species_count_dict = {}
species_list = []

for file in var_level_files:
    ht = hl.read_table(file)

    species = file.split('/')[-1].removesuffix('_lifted_vars.ht')
    species_list.append(species)
    call_count = str([dtype for name, dtype in ht.row.dtype.items() if name == 'call_data'][0]).count('call')
    species_count_dict[species] = call_count
    
    all_positions = all_positions.annotate(
        species = ht[all_positions.key][f'{species}_data']
    )
    all_positions = all_positions.rename({'species': f'{species}_data'})
    
    all_positions = all_positions.checkpoint(f'./tmp/{species}_data_2.ht', overwrite=True)

species_list = list(species_count_dict.keys())

species_count_dict = hl.literal(species_count_dict)

# Compute species-level summaries at the row level
all_positions = all_positions.annotate(
    species_any_locus_defined=hl.sum([
        all_positions[f'{species}_data'].any(lambda x: hl.is_defined(x.species_locus)) for species in species_list
    ]),
    species_any_alleles_defined=hl.sum([
        all_positions[f'{species}_data'].any(lambda x: hl.is_defined(x.species_alleles)) for species in species_list
    ])
)

all_positions = all_positions.checkpoint('all_positions_cpt2.ht', overwrite=True)

all_positions = all_positions.annotate(
    AC_sum_mean=hl.sum(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.AC)) for species in species_list])),
    AN_sum_mean=hl.sum(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.AN)) for species in species_list])),
    AC_sum_min=hl.sum(hl.array([hl.min(all_positions[f'{species}_data'].map(lambda x: x.AC)) for species in species_list])),
    AN_sum_for_min=hl.sum(hl.array([all_positions[f'{species}_data'][hl.argmin(all_positions[f'{species}_data'].map(lambda x: x.AC))].AN for species in species_list])),
    AC_sum_max=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.AC)) for species in species_list])),
    AN_sum_for_max=hl.sum(hl.array([all_positions[f'{species}_data'][hl.argmax(all_positions[f'{species}_data'].map(lambda x: x.AC))].AN for species in species_list]))
)

all_positions = all_positions.checkpoint('all_positions_cpt3.ht', overwrite=True)

all_positions = all_positions.annotate(
    species=hl.array(hl.set(hl.flatten(hl.set([hl.flatten(hl.set(all_positions[f'{species}_data'].species)) for species in species_list])))),
    family=hl.array(hl.set(hl.flatten(hl.set([hl.flatten(hl.set(all_positions[f'{species}_data'].family)) for species in species_list])))),
    reference=hl.array(hl.set(hl.flatten(hl.set([hl.flatten(hl.set(all_positions[f'{species}_data'].reference)) for species in species_list])))),
    order=hl.array(hl.set(hl.flatten(hl.set([hl.flatten(hl.set(all_positions[f'{species}_data'].order)) for species in species_list])))),
    n_hom_sum_max=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.n_hom)) for species in species_list])),
    n_hom_sum_min=hl.sum(hl.array([hl.min(all_positions[f'{species}_data'].map(lambda x: x.n_hom)) for species in species_list])),
    n_hom_sum_mean=hl.sum(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.n_hom)) for species in species_list]))
)

all_positions = all_positions.checkpoint('all_positions_cpt4.ht', overwrite=True)

all_positions = all_positions.annotate(
    mean_min_minhamming=hl.mean(hl.array([hl.min(all_positions[f'{species}_data'].map(lambda x: x.min_hamming)) for species in species_list])),
    mean_max_minhamming=hl.mean(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.min_hamming)) for species in species_list])),
    mean_pentamer_error=hl.mean(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.pentamer_error)) for species in species_list])),
    mean_max_pentamer_error=hl.mean(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.pentamer_error)) for species in species_list])),
    mean_min_pentamer_error=hl.mean(hl.array([hl.min(all_positions[f'{species}_data'].map(lambda x: x.pentamer_error)) for species in species_list])),
    mean_human2species_mappings=hl.mean(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.human2species_mappings)) for species in species_list])),
    mean_overlapping_annotations=hl.mean(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.overlapping_annotations)) for species in species_list]))
)

all_positions = all_positions.checkpoint('all_positions_cpt5.ht', overwrite=True)

all_positions = all_positions.annotate(
    sum_max_Platyrrhini_hom=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Platyrrhini_hom)) for species in species_list])),
    sum_max_Catarrhini_hom=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Catarrhini_hom)) for species in species_list])),
    sum_max_Strepsirrhini_hom=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Strepsirrhini_hom)) for species in species_list])),
    sum_max_Tarsiiformes_hom=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Tarsiiformes_hom)) for species in species_list])),
    sum_max_Platyrrhini_het=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Platyrrhini_het)) for species in species_list])),
    sum_max_Catarrhini_het=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Catarrhini_het)) for species in species_list])),
    sum_max_Strepsirrhini_het=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Strepsirrhini_het)) for species in species_list])),
    sum_max_Tarsiiformes_het=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: x.Tarsiiformes_het)) for species in species_list]))
)

all_positions = all_positions.checkpoint('all_positions_cpt6.ht', overwrite=True)

all_positions = all_positions.annotate(
    mean_ref_match=hl.mean(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.ref_match)) for species in species_list])),
    sum_max_hom_ref_match=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: hl.if_else((x.ref_match == 1), x.n_hom, 0))) for species in species_list])),
    sum_max_het_ref_match=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: hl.if_else((x.ref_match == 1), x.n_het, 0))) for species in species_list])),
    sum_max_hom_ref_mismatch=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: hl.if_else((x.ref_match == 0), x.n_hom, 0))) for species in species_list])),
    sum_max_het_ref_mismatch=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: hl.if_else((x.ref_match == 0), x.n_het, 0))) for species in species_list])),
    sum_max_AN_ref_match=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: hl.if_else((x.ref_match == 1), x.AN, 0))) for species in species_list])),
    sum_max_AN_ref_mismatch=hl.sum(hl.array([hl.max(all_positions[f'{species}_data'].map(lambda x: hl.if_else((x.ref_match == 0), x.AN, 0))) for species in species_list]))
)

all_positions = all_positions.checkpoint('all_positions_cpt7.ht', overwrite=True)

all_positions = all_positions.annotate(
    AF_mean=hl.if_else(all_positions.AN_sum_mean > 0, all_positions.AC_sum_mean/all_positions.AN_sum_mean, 0.0),
    AF_min=hl.if_else(all_positions.AN_sum_for_min > 0, all_positions.AC_sum_min/all_positions.AN_sum_for_min, 0.0),
    AF_max=hl.if_else(all_positions.AN_sum_for_max > 0, all_positions.AC_sum_max/all_positions.AN_sum_for_max, 0.0),
    mean_human2species_mappings=hl.mean(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.human2species_mappings)) for species in species_list])),
    mean_overlapping_annotations=hl.mean(hl.array([hl.mean(all_positions[f'{species}_data'].map(lambda x: x.overlapping_annotations)) for species in species_list])),
)

all_positions = all_positions.annotate(
    weighted_pentamer_error = hl.array([hl.if_else(
        # Condition: Check if the sum of AN for valid entries is greater than 0
        hl.sum(
            all_positions[f'{species}_data']
            .filter(lambda x: hl.is_defined(x.pentamer_error) & hl.is_defined(x.AN))
            .map(lambda x: x.AN)
        ) > 0,
        # True branch: Calculate the weighted average
        hl.sum(
            all_positions[f'{species}_data']
            .filter(lambda x: hl.is_defined(x.pentamer_error) & hl.is_defined(x.AN) )
            .map(lambda x: x.pentamer_error * x.AN)
        ) / hl.sum(
            all_positions[f'{species}_data']
            .filter(lambda x: hl.is_defined(x.pentamer_error) & hl.is_defined(x.AN))
            .map(lambda x: x.AN)
        ) * species_count_dict[species],
        # False branch: Return missing if sum_an is 0 or missing
            hl.missing(hl.tfloat64)
        )
        for species in species_list
    ])
)

all_positions = all_positions.annotate(
    weighted_mean_pentamer_error = hl.sum(all_positions.weighted_pentamer_error) / hl.sum(species_count_dict.values())
)

all_positions = all_positions.checkpoint('aggregated_primate_counts.ht', overwrite=True)