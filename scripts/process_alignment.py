import sys
import re
import gzip
import glob
from Bio import AlignIO, SeqIO
import polars as pl
from joblib import Parallel, delayed
import random
import os

def make_alias_dict(alias_file_path):
    with open(alias_file_path, 'r') as file:
        lines = file.readlines()

    headers = lines[0].replace('# ', '').strip().split('\t')
    return {header: [row.split('\t')[i] for row in lines[1:]] for i, header in enumerate(headers)}

def determine_name_format(file_path, alias_dict, species=None, file_type='maf'):
    if file_type == 'maf':
        alignments = AlignIO.parse(file_path, 'maf')
        sequence_name = None
        for alignment in alignments:
            for record in alignment:
                if record.id.startswith(species):
                    parts = record.id.split('.')
                    if len(parts) > 1:
                        sequence_name = parts[1]
                        if sequence_name == 'ContigUN':
                            continue
                        elif sequence_name:
                            break
        if not sequence_name:
            return None
    else:  # ref
        with gzip.open(file_path, 'rt') as fna_file:
            sequence_name = next(SeqIO.parse(fna_file, 'fasta')).description.split(' ')[0].split('.')[0]

    for key, value in alias_dict.items():
        if sequence_name in [item.split('.')[0] for item in value]:
            return key
    return None

def map_aliases(alias_file_path, maf_naming, ref_naming):
    with open(alias_file_path, 'r') as file:
        lines = file.readlines()
    headers = lines[0].replace('# ', '').strip().split('\t')
    maf_column, ref_column = headers.index(maf_naming), headers.index(ref_naming)

    mapping = {}
    for line in lines[1:]:
        ref, maf = line.strip().split('\t')[ref_column], line.strip().split('\t')[maf_column]
        ref_short, maf_short = ref.split('.')[0], maf.split('.')[0]
        mapping.update({
            ref: maf,
            ref_short: maf
        })
    return mapping

def load_contig_ids(ref_file_path, maf_filepath, alias_file_path, species):
    alias_dict = make_alias_dict(alias_file_path)
    maf_naming = determine_name_format(maf_filepath, alias_dict, species, 'maf')
    ref_naming = determine_name_format(ref_file_path, alias_dict, file_type='ref')

    print(f"Naming Formats:\n{maf_naming}\n{ref_naming}")

    if not (maf_naming and ref_naming):
        raise ValueError(f"Naming Mismatch:\n{maf_filepath}\n{ref_file_path}")

    contig_ids = {}
    
    alias_mapping = (lambda x: x) if maf_naming == ref_naming else map_aliases(alias_file_path, maf_naming, ref_naming).get

    with gzip.open(ref_file_path, 'rt') as fna_file:
        for record in SeqIO.parse(fna_file, 'fasta'):
            contigs = []
            contig = re.search(r"\bContig\d+\b", record.description)
            if contig:
                contigs.append(contig.group())
            contigs.append(record.description.split(' ')[0].split('.')[0])
            for contig in contigs:
                key = alias_mapping(contig)
                if key:
                    contig_ids[key] = record.seq
                    contig_ids[key.split('.')[0]] = record.seq    

    return contig_ids

def get_alignment_orig(alignment):
    length = alignment.get_alignment_length()

    hg38_record = None
    nhp_record = None

    for record in alignment:
        if record.id.startswith('hg38'):
            hg38_record = record
        elif record.id.startswith(species):
            nhp_record = record
    if nhp_record is None:
        return None
    else:
        alignment_df = pl.DataFrame({
            'hg38_chr': [hg38_record.id.split('.')[1]] * length,
            'hg38_start': [hg38_record.annotations['start']] * length,
            'hg38_seq': list(hg38_record.seq.upper()),
            'hg38_strand': [hg38_record.annotations['strand']] * length,
            'hg38_srcsize': [hg38_record.annotations['srcSize']] * length,
            f'{species}_contig': [nhp_record.id.split('.')[1]] * length,
            f'{species}_start': [nhp_record.annotations['start']] * length,
            f'{species}_seq_msa': list(nhp_record.seq.upper()),
            f'{species}_strand': [nhp_record.annotations['strand']] * length,
            f'{species}_srcsize': [nhp_record.annotations['srcSize']] * length,
            'percent_aligned': [nhp_record.annotations['size']/hg38_record.annotations['size']] * length,
        })

    return alignment_df.with_columns([
        pl.when(pl.col('hg38_seq') != "-").then(1).otherwise(0).alias('hg38_nongap'),
        pl.when(pl.col(f'{species}_seq_msa') != "-").then(1).otherwise(0).alias(f'{species}_nongap')
    ]).with_columns([
        pl.col('hg38_nongap').cum_sum().alias('hg38_cumsum'),
        pl.col(f'{species}_nongap').cum_sum().alias(f'{species}_cumsum'),
    ])


def get_ref_base(pos, contig, contig_dict):
    if pos == -1:
        return "-"
    else:
        try:
            result = contig_dict[contig][pos].upper()
            return result
        except KeyError as e:
            if contig not in contig_error_list:
                contig_error_list.append(contig)
            return "-"

def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(dna_sequence.upper()))

if __name__ == "__main__":
    # Get command-line arguments
    gca_id, species, maf_file_path, alias_file_path, ref_file_path, output_dir = sys.argv[1:7]

    print(f"Processing maf for {species}")

    id = gca_id.split('_')[1].split('.')[0]
    
    maf_basename = os.path.splitext(os.path.basename(maf_file_path))[0]

    Contig_IDs = load_contig_ids(ref_file_path, maf_file_path, alias_file_path, species)

    alignments = AlignIO.parse(maf_file_path, 'maf')

    results = Parallel(n_jobs=-1)(delayed(get_alignment_orig)(alignment) for alignment in alignments)
    results_len = len(results)
    results = [result for result in results if result is not None]
    #If result is None it means there was no NHP record for that alignment (i.e. the NHP did not align)
    print(f'{results_len - len(results)} of {results_len} alignments did not align to hg38')

    all_alignments = pl.concat(results, how="vertical_relaxed")

    all_alignments = (all_alignments
        .lazy()
        .with_columns([
            pl.when(pl.col("hg38_nongap") == 1)
              .then(pl.col("hg38_cumsum") + pl.col("hg38_start") - 1)
              .otherwise(None)
              .alias("hg38_position"),
            pl.when((pl.col(f'{species}_nongap') == 1) & (pl.col(f'{species}_strand') == -1))
              .then(pl.col(f'{species}_srcsize') - pl.col(f'{species}_start') - pl.col(f'{species}_cumsum'))
              .otherwise(
                pl.when((pl.col(f'{species}_nongap') == 1) & (pl.col(f'{species}_strand') == 1))
                  .then(pl.col(f'{species}_srcsize') - pl.col(f'{species}_start') + pl.col(f'{species}_cumsum'))
              )
              .alias(f'{species}_position_rev'),
            pl.when((pl.col(f'{species}_nongap') == 1) & (pl.col(f'{species}_strand') == -1))
              .then(pl.col(f'{species}_start') + pl.col(f'{species}_cumsum') * -1)
              .otherwise(
                pl.when((pl.col(f'{species}_nongap') == 1) & (pl.col(f'{species}_strand') == 1))
                  .then(pl.col(f'{species}_start') + pl.col(f'{species}_cumsum') - 1)
              )
              .alias(f'{species}_position_orig')]).with_columns(
            pl.when(pl.col(f'{species}_nongap') == 1)
                .then(
                    pl.when(pl.col(f'{species}_strand') == 1)
                    .then(pl.col(f'{species}_position_orig'))
                    .when(pl.col(f'{species}_strand') == -1)
                    .then(pl.col(f'{species}_position_rev'))
                    .otherwise(-1))
                .otherwise(-1)
                .alias(f'{species}_position')
        )
        .filter(pl.col("hg38_position").is_not_null())
        .drop(["hg38_start", "hg38_nongap", "hg38_cumsum", "hg38_srcsize",
               f'{species}_start', f'{species}_cumsum', f'{species}_srcsize'])
        .collect()
    )
    
    contig_error_list = []

    ref_bases = all_alignments.map_rows(lambda row: get_ref_base(
            row[11],
            row[3],
            Contig_IDs
        )).to_series()
    
    print(contig_error_list)

    all_alignments = all_alignments.with_columns([
            ref_bases.alias('seq_refs')
        ]).with_columns([
            pl.when((pl.col(f'{species}_strand') == -1) & (pl.col(f'{species}_nongap') == 1))
              .then(pl.col('seq_refs').map_elements(lambda x: reverse_complement(x), return_dtype=pl.Utf8))
              .otherwise(pl.col('seq_refs'))
              .alias('aligned_base')
        ]).with_columns([
                pl.when(pl.col('aligned_base') == pl.col(f'{species}_seq_msa'))
                  .then(1)
                  .otherwise(0)
                  .alias('match')
        ])

    all_alignments.write_parquet(f'{output_dir}/{maf_basename}_alignment.parquet')
