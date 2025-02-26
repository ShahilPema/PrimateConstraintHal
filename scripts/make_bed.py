#!usr/bin/env python

import polars as pl
import re
import argparse

def validate_and_fix_contig(contig):
    contig = re.sub(r'^chr', '', contig, flags=re.IGNORECASE)
    if re.match(r'^((?:[1-9]|1[0-9]|2[0-2])|X|Y)$', contig, re.IGNORECASE):
        return f"chr{contig.upper()}"
    return None

def parse_gff(mane_transcript_path, gff_path):
    with open(mane_transcript_path, 'r') as file:
        mane_transcripts = set(line.strip() for line in file)

    elements = []
    lnc_RNA = None

    with open(gff_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                fields = line.split('\t')
                if fields[2] in ['lnc_RNA', 'miRNA', 'scRNA', 'snRNA', 'snoRNA', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'exon']:
                    chrom = validate_and_fix_contig(fields[0])
                    if chrom is None:
                        continue

                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    region_type = 'CDS' if fields[2] == 'CDS' else '5PrimeUTR' if fields[2] == 'five_prime_UTR' else '3PrimeUTR' if fields[2] == 'three_prime_UTR' else fields[2]

                    transcript_id = next(item.split(':')[1] for item in fields[8].split(';') if 'transcript:' in item).strip().split('.')[0]

                    if region_type in ['CDS', '5PrimeUTR', '3PrimeUTR'] and transcript_id not in mane_transcripts:
                        continue

                    if region_type == 'lnc_RNA':
                        if any('Ensembl_canonical' in item for item in fields[8].split(';')):
                            lnc_RNA = transcript_id
                            continue
                        else:
                            continue

                    if region_type == 'exon':
                        if transcript_id == lnc_RNA:
                            elements.append([chrom, start - 1, end, strand, f"{transcript_id}_lncRNA"])
                            if any(item == 'rank=1' for item in fields[8].split(';')):
                                if strand == '+':
                                    elements.append([chrom, start - 501, start - 1, strand, f"{transcript_id}_lncRNAPromoterProximal500"])
                                    elements.append([chrom, start - 1001, start - 501, strand, f"{transcript_id}_lncRNAPromoterDistal500"])
                                else:
                                    elements.append([chrom, end, end + 500, strand, f"{transcript_id}_lncRNAPromoterProximal500"])
                                    elements.append([chrom, end + 500, end + 1000, strand, f"{transcript_id}_lncRNAPromoterDistal500"])
                            continue
                        else:
                            continue
                    else:
                        elements.append([chrom, start - 1, end, strand, f"{transcript_id}_{region_type}"])

    return elements

def process_5prime_utrs(data, strand):
    data = data.with_columns((pl.col('end') - pl.col('start')).alias('length'))
    dataframes = data.partition_by('region')
    processed_dataframes = []

    for df in dataframes:
        if strand == '+':
            df = df.sort('start')
        else:
            df = df.sort('start', descending=True)        
        df = df.with_columns(
            pl.col('length').cum_sum().alias('cum_sum')).with_columns(
            (pl.lit(50) - pl.col('cum_sum')).alias('remaining')
        )

        if strand == '+':
            df = df.with_columns(
                pl.when(pl.col('cum_sum') > 50)
                .then(pl.col('start') + pl.col('remaining') + pl.col('length'))
                .otherwise(pl.col('end'))
                .alias('end')
            )
        else:
            df = df.with_columns(
                pl.when(pl.col('cum_sum') > 50)
                .then(pl.col('end') - (pl.col('remaining') + pl.col('length')))
                .otherwise(pl.col('start'))
                .alias('start')
            )

        df = df.filter((pl.col('remaining') + pl.col('length')) > 0).with_columns(
            (pl.col('end') - pl.col('start')).alias('length'),
            pl.col('region').str.split('_').list.get(0).alias('transcript'),
            pl.lit('5primeUTR50bp').alias('region_type')
        ).with_columns(
            (pl.col('transcript') + pl.lit('_') + pl.col('region_type')).alias('region')
        )

        processed_dataframes.append(df)

    return pl.concat(processed_dataframes)

def generate_promoter_data(utrs_50bp, strand):
    if strand == '+':
        utrs_50bp = utrs_50bp.sort('start')
    else:
        utrs_50bp = utrs_50bp.sort('start', descending=True)       
    promoter_data = utrs_50bp.group_by('transcript').first()

    if strand == '+':
        promoter_proximal = promoter_data.with_columns(
            pl.col('start').alias('end'),
            (pl.col('start') - pl.lit(500)).alias('start'),
            pl.lit('PromoterProximal500').alias('region_type')
        )
        promoter_distal = promoter_data.with_columns(
            (pl.col('start') - pl.lit(500)).alias('end'),
            (pl.col('start') - pl.lit(1000)).alias('start'),
            pl.lit('PromoterDistal500').alias('region_type')
        )
    else:
        promoter_proximal = promoter_data.with_columns(
            pl.col('end').alias('start'),
            (pl.col('end') + pl.lit(500)).alias('end'),
            pl.lit('PromoterProximal500').alias('region_type')
        )
        promoter_distal = promoter_data.with_columns(
            (pl.col('end') + pl.lit(500)).alias('start'),
            (pl.col('end') + pl.lit(1000)).alias('end'),
            pl.lit('PromoterDistal500').alias('region_type')
        )

    promoter_proximal = promoter_proximal.with_columns(
        (pl.col('end') - pl.col('start')).alias('length'),
        (pl.col('transcript') + pl.lit('_') + pl.col('region_type')).alias('region')
    ).select(['hg38_chr', 'start', 'end', 'strand', 'region'])

    promoter_distal = promoter_distal.with_columns(
        (pl.col('end') - pl.col('start')).alias('length'),
        (pl.col('transcript') + pl.lit('_') + pl.col('region_type')).alias('region')
    ).select(['hg38_chr', 'start', 'end', 'strand', 'region'])

    return promoter_proximal, promoter_distal

def main(mane_transcript_path, gff_path, output_path):
    regions = parse_gff(mane_transcript_path, gff_path)
    bed_data = pl.DataFrame(regions, schema={'hg38_chr': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64, 'strand': pl.Utf8, 'region': pl.Utf8}, orient='row')
    bed_data.write_csv('features_by_loci.bed', separator='\t')

    data = pl.read_csv('features_by_loci.bed', separator='\t')
    pos_strand_5prime = data.filter((pl.col('strand') == "+") & (pl.col('region').str.contains('5PrimeUTR')))
    neg_strand_5prime = data.filter((pl.col('strand') == "-") & (pl.col('region').str.contains('5PrimeUTR')))

    pos_utrs_50bp = process_5prime_utrs(pos_strand_5prime, '+')
    neg_utrs_50bp = process_5prime_utrs(neg_strand_5prime, '-')

    pos_promoter_proximal, pos_promoter_distal = generate_promoter_data(pos_utrs_50bp, '+')
    neg_promoter_proximal, neg_promoter_distal = generate_promoter_data(neg_utrs_50bp, '-')

    pos_utrs_50bp = pos_utrs_50bp.select(['hg38_chr', 'start', 'end', 'strand', 'region'])
    neg_utrs_50bp = neg_utrs_50bp.select(['hg38_chr', 'start', 'end', 'strand', 'region'])

    processed_dataframe = pl.concat([
        data,
        pos_utrs_50bp, pos_promoter_proximal, pos_promoter_distal,
        neg_utrs_50bp, neg_promoter_proximal, neg_promoter_distal
    ])

    processed_dataframe.write_csv(output_path, separator='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process GFF file and create BED file')
    parser.add_argument('--mane_transcript', required=True, help='Path to MANE transcript file')
    parser.add_argument('--gff_file', required=True, help='Path to GFF file')
    parser.add_argument('--output', required=True, help='Output BED file path')
    
    args = parser.parse_args()
    
    main(args.mane_transcript, args.gff_file, args.output)
