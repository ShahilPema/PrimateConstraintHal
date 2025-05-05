#!/usr/bin/env python3

import hail as hl
import os
import polars as pl
import sys
import psutil
import argparse
import math
from typing import Optional

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process a VCF file with Hail.")
    parser.add_argument("--cpus", required=True, help="Total cpus on machine")
    parser.add_argument("--memory", required=True, help="Total memory on machine")
    parser.add_argument("--species", required=True, help="Species name (e.g., 'Macaca_mulatta').")
    parser.add_argument("--vcf_file", required=True, help="Path to the VCF file.")
    parser.add_argument("--alias_file", required=True, help="Path to the alias file.")
    parser.add_argument("--vcf_format", required=True, help="VCF contig format (e.g., 'ucsc').")
    parser.add_argument("--hal_format", required=True, help="HAL contig format (e.g., 'refseq').")
    parser.add_argument("--primate_fasta_path", required=True, help="Path to the FASTA file for the primate genome.")
    parser.add_argument("--primate_index_path", required=True, help="Path to the FASTA index file (.fai).")
    parser.add_argument("--tmpdir", required=True, help="Temporary directory")
    parser.add_argument("--output", required=True, help="Output directory")
    args = parser.parse_args()

    config = {
        'spark.driver.memory': f'{args.memory}',  #Set to total memory
        'spark.executor.memory': f'{args.memory}',
        'spark.local.dir': args.tmpdir,
        'spark.driver.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}',
        'spark.executor.extraJavaOptions': f'-Djava.io.tmpdir={args.tmpdir}'
    }
 
    hl.init(spark_conf=config, master=f'local[{args.cpus}]', tmp_dir=args.tmpdir, local_tmpdir=args.tmpdir)

    # Load alias file and create mapping
    if args.vcf_format != args.hal_format:
        alias_data = pl.read_csv(args.alias_file, separator='\t')
        alias_data.columns = [item.replace(" ", "").replace("#", "") for item in alias_data.columns]
        mapping_data = alias_data.select(args.vcf_format, args.hal_format)
        mapping_dict = dict(mapping_data.iter_rows())
        mapping_dict = {k: v for k, v in mapping_dict.items() if k is not None and v is not None}
    else:
        mapping_dict=None

    # Create reference genome
    primate_rg = hl.ReferenceGenome.from_fasta_file(
        name=args.species,
        fasta_file=args.primate_fasta_path,
        index_file=args.primate_index_path
    )

    # Process VCF file
    def process_vcf(vcf_file: str, species: str, mapping_dict: Optional[dict], cpus: int) -> hl.MatrixTable:
        mt = hl.import_vcf(
            vcf_file,
            reference_genome=species,
            contig_recoding=mapping_dict,
            force_bgz=True,
            min_partitions=cpus
        )
        
        mt = mt.annotate_rows(
            filter_stats=hl.struct(
                fail_gatk_filters=(
                    (mt.info.QD < 2.0) |
                    (mt.info.FS > 60.0) |
                    (mt.info.MQ < 40.0) |
                    (mt.info.SOR > 3.0) |
                    (hl.is_defined(mt.info.ReadPosRankSum) & (mt.info.ReadPosRankSum < -8.0)) |
                    (hl.is_defined(mt.info.MQRankSum) & (mt.info.MQRankSum < -12.5))
                ),
                fail_het_depth=hl.agg.any(
                    mt.GT.is_het() & ((mt.AD[0] < 3) | (mt.AD[1] < 3))
                )
            )
        )
        mt = mt.annotate_cols(
            pass_filter=hl.agg.count_where(~mt.filter_stats.fail_gatk_filters & ~mt.filter_stats.fail_het_depth),
            total=hl.agg.count(),
            het_var=hl.agg.count_where(mt.GT.is_het_non_ref()),
            het_ref=hl.agg.count_where(mt.GT.is_het_ref()),
            hom_var=hl.agg.count_where(mt.GT.is_hom_var())
        )
        
        mt = mt.select_rows()
        
        mt = mt.select_entries(mt.GT)
        return mt

    basename = os.path.basename(args.vcf_file).split('.')[0]
    mt = process_vcf(args.vcf_file, args.species, mapping_dict, 16)
    mt.write(f'{args.output}/{basename}.mt', overwrite=True)

if __name__ == "__main__":
    main()

