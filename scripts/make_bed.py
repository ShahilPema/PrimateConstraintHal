#!/usr/bin/env python

import polars as pl
import re
import argparse
import os

def validate_and_fix_contig(contig):
    contig = re.sub(r'^chr', '', contig, flags=re.IGNORECASE)
    if re.match(r'^((?:[1-9]|1[0-9]|2[0-2])|X|Y)$', contig, re.IGNORECASE):
        return f"chr{contig.upper()}"
    return None

def parse_gff(gff_path):
    elements = []
    selected_lnc_RNAs = set()
    
    # Dictionary to collect exons by transcript for splice site identification
    transcript_exons = {}
    # Dictionary to track lncRNA transcripts by gene
    lnc_RNA_by_gene = {}
    transcript_exon_lengths = {}
    
    # Dictionaries to store UTR and CDS regions for each transcript
    five_prime_utrs = {}
    three_prime_utrs = {}
    cds_regions = {}
    
    # Track processed transcripts to avoid duplicate processing
    processed_transcripts = set()
    # Track potential lncRNA transcripts
    potential_lnc_RNAs = set()
    # Track MANE Select transcripts
    mane_transcripts = set()
    
    # Store lncRNA transcript info for tag tracking
    lnc_RNA_transcript_info = {}

    # First perform a complete scan of the GFF file to collect all relevant information
    with open(gff_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                fields = line.split('\t')
                feature_type = fields[2]
                
                # Skip invalid contigs early
                chrom = validate_and_fix_contig(fields[0])
                if chrom is None:
                    continue
                
                # Extract transcript ID
                transcript_id = next((item.split(':')[1] for item in fields[8].split(';') 
                                   if 'transcript:' in item), '').strip().split('.')[0]
                
                if not transcript_id:
                    continue

                # Check if this is a MANE Select transcript
                if feature_type == 'mRNA':
                    tag_field = next((field for field in fields[8].split(';') if field.startswith('tag=')), '')
                    if 'MANE_Select' in tag_field:
                        mane_transcripts.add(transcript_id)
                
                # Check if this is a MANE transcript
                is_mane = transcript_id in mane_transcripts
                
                # Process lncRNA feature types
                if feature_type in ['lnc_RNA', 'lncRNA', 'lincRNA', 'long_non_coding_RNA']:
                    potential_lnc_RNAs.add(transcript_id)
                    
                    # For lncRNAs, also extract gene ID and tags for later selection
                    gene_id = next((item.split(':')[1] for item in fields[8].split(';') 
                                  if 'Parent=gene:' in item), '').strip().split('.')[0]
                    
                    if gene_id:
                        tag_field = next((field for field in fields[8].split(';') if field.startswith('tag=')), '')
                        tags = tag_field.split('=')[1].split(',') if tag_field else []
                        is_gencode_basic = 'gencode_basic' in tags
                        is_ensembl_canonical = 'Ensembl_canonical' in tags
                        
                        # Store transcript info for later selection
                        lnc_RNA_transcript_info[transcript_id] = (gene_id, is_gencode_basic, is_ensembl_canonical)
                
                # Process exons for all transcripts
                if feature_type == 'exon':
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    exon_length = end - start + 1
                    
                    # Store exon length for all transcripts (needed for lncRNA selection)
                    if transcript_id not in transcript_exon_lengths:
                        transcript_exon_lengths[transcript_id] = []
                    transcript_exon_lengths[transcript_id].append(exon_length)
                    
                    # Store full exon info for MANE and potential lncRNAs
                    if is_mane or transcript_id in potential_lnc_RNAs:
                        # Get exon rank if available
                        exon_rank = next((int(item.split('=')[1]) for item in fields[8].split(';') 
                                        if item.startswith('rank=')), 0)
                        
                        if transcript_id not in transcript_exons:
                            transcript_exons[transcript_id] = []
                            
                        transcript_exons[transcript_id].append({
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'rank': exon_rank
                        })
                
                # Process UTR and CDS regions (for MANE transcripts only)
                if is_mane:
                    if feature_type == 'five_prime_UTR':
                        if transcript_id not in five_prime_utrs:
                            five_prime_utrs[transcript_id] = []
                        five_prime_utrs[transcript_id].append({'start': int(fields[3]), 'end': int(fields[4])})
                    
                    elif feature_type == 'three_prime_UTR':
                        if transcript_id not in three_prime_utrs:
                            three_prime_utrs[transcript_id] = []
                        three_prime_utrs[transcript_id].append({'start': int(fields[3]), 'end': int(fields[4])})
                    
                    elif feature_type == 'CDS':
                        if transcript_id not in cds_regions:
                            cds_regions[transcript_id] = []
                        cds_regions[transcript_id].append({'start': int(fields[3]), 'end': int(fields[4])})

    # Organize lncRNA transcripts by gene
    for transcript_id, (gene_id, is_basic, is_canonical) in lnc_RNA_transcript_info.items():
        if transcript_id in transcript_exon_lengths and is_basic:
            total_exon_length = sum(transcript_exon_lengths[transcript_id])
            
            if gene_id not in lnc_RNA_by_gene:
                lnc_RNA_by_gene[gene_id] = []
                
            lnc_RNA_by_gene[gene_id].append((transcript_id, is_basic, is_canonical, total_exon_length))
    
    # Now select transcripts based on complete gene information
    for gene_id, transcripts in lnc_RNA_by_gene.items():
        # First look for canonical+basic transcript
        canonical_basic = [t for t in transcripts if t[1] and t[2]]  # is_basic and is_canonical
        if canonical_basic:
            # Add the canonical+basic transcript
            selected_lnc_RNAs.add(canonical_basic[0][0])
        else:
            # If no canonical+basic, get the longest basic transcript
            basic_transcripts = [t for t in transcripts if t[1]]  # is_basic
            if basic_transcripts:
                longest_basic = max(basic_transcripts, key=lambda x: x[3])  # sort by length
                selected_lnc_RNAs.add(longest_basic[0])

    # Remove exon information for non-selected lncRNAs to save memory
    for transcript_id in list(transcript_exons.keys()):
        if transcript_id in potential_lnc_RNAs and transcript_id not in selected_lnc_RNAs:
            del transcript_exons[transcript_id]

    # Process collected exons for splice sites
    for transcript_id, exons in transcript_exons.items():
        # Skip if we've already processed this transcript
        if transcript_id in processed_transcripts:
            continue
            
        processed_transcripts.add(transcript_id)
        
        # Sort exons by rank or position if rank not available
        if all(exon['rank'] for exon in exons):
            exons.sort(key=lambda x: x['rank'])
        else:
            exons.sort(key=lambda x: x['start'])
            
        # Remove duplicate exons (same start and end)
        unique_exons = []
        seen_positions = set()
        
        for exon in exons:
            pos_key = (exon['start'], exon['end'])
            if pos_key not in seen_positions:
                seen_positions.add(pos_key)
                unique_exons.append(exon)
                
        exons = unique_exons
        num_exons = len(exons)
        
        if num_exons <= 1:  # Skip single-exon transcripts
            continue
            
        # Skip transcript if any exons overlap (strand-aware check)
        has_invalid_exons = False
        strand = exons[0]['strand']
        for i in range(num_exons - 1):
            # Different overlap check based on strand
            if strand == '-':
                # For minus strand: if next_end >= current_start, they overlap
                if exons[i+1]['end'] >= exons[i]['start']:
                    has_invalid_exons = True
                    break
            else:
                # For plus strand: if current_end >= next_start, they overlap
                if exons[i]['end'] >= exons[i+1]['start']:
                    has_invalid_exons = True
                    break
                    
        if has_invalid_exons:
            continue
            
        # Get transcript type
        is_lnc_RNA = transcript_id in selected_lnc_RNAs
        
        # Get all splice sites for this transcript
        total_exons = len(exons)
        chrom = exons[0]['chrom']
        strand = exons[0]['strand']
        
        # Process splice sites differently based on strand
        if strand == '-':
            # For minus strand, process exons in descending order
            for i in range(total_exons - 1):
                current_exon = exons[i]
                next_exon = exons[i+1]
                
                # For minus strand: if next_end >= current_start, they overlap
                if next_exon['end'] >= current_exon['start']:
                    continue
                
                # Donor site (end of current exon on minus strand)
                donor_pos = current_exon['end']
                splice_site_type = "SpliceSite"  # Default label
                
                if not is_lnc_RNA:  # Only check UTRs for protein-coding genes
                    # Check if donor site is in 3' UTR (remember strand inversion)
                    if transcript_id in three_prime_utrs:
                        for utr in three_prime_utrs[transcript_id]:
                            if donor_pos >= utr['start'] and donor_pos <= utr['end']:
                                splice_site_type = "3PrimeUTRSpliceSite"
                                break
                    
                    # Check if donor site is in 5' UTR if not found in 3' UTR
                    if splice_site_type == "SpliceSite" and transcript_id in five_prime_utrs:
                        for utr in five_prime_utrs[transcript_id]:
                            if donor_pos >= utr['start'] and donor_pos <= utr['end']:
                                splice_site_type = "5PrimeUTRSpliceSite"
                                break
                
                # Add donor splice site with appropriate label
                if is_lnc_RNA:
                    elements.append([chrom, donor_pos, donor_pos + 2, strand, f"{transcript_id}-lncRNAspliceSite"])
                else:
                    elements.append([chrom, donor_pos, donor_pos + 2, strand, f"{transcript_id}-{splice_site_type}"])
                
                # Acceptor site (start of next exon on minus strand)
                acceptor_pos = next_exon['start']
                splice_site_type = "SpliceSite"  # Reset label
                
                if not is_lnc_RNA:  # Only check UTRs for protein-coding genes
                    # Check if acceptor site is in 3' UTR (remember strand inversion)
                    if transcript_id in three_prime_utrs:
                        for utr in three_prime_utrs[transcript_id]:
                            if acceptor_pos >= utr['start'] and acceptor_pos <= utr['end']:
                                splice_site_type = "3PrimeUTRSpliceSite"
                                break
                    
                    # Check if acceptor site is in 5' UTR if not found in 3' UTR
                    if splice_site_type == "SpliceSite" and transcript_id in five_prime_utrs:
                        for utr in five_prime_utrs[transcript_id]:
                            if acceptor_pos >= utr['start'] and acceptor_pos <= utr['end']:
                                splice_site_type = "5PrimeUTRSpliceSite"
                                break
                
                # Add acceptor splice site with appropriate label
                if is_lnc_RNA:
                    elements.append([chrom, acceptor_pos - 3, acceptor_pos - 1, strand, f"{transcript_id}-lncRNAspliceSite"])
                else:
                    elements.append([chrom, acceptor_pos - 3, acceptor_pos - 1, strand, f"{transcript_id}-{splice_site_type}"])
                
        else:  # plus strand
            # For plus strand, process exons in ascending order
            for i in range(total_exons - 1):
                current_exon = exons[i]
                next_exon = exons[i+1]
                
                # For plus strand: if current_end >= next_start, they overlap
                if current_exon['end'] >= next_exon['start']:
                    continue
                
                # Donor site (end of current exon)
                donor_pos = current_exon['end']
                splice_site_type = "SpliceSite"  # Default label
                
                if not is_lnc_RNA:  # Only check UTRs for protein-coding genes
                    # Check if donor site is in 5' UTR
                    if transcript_id in five_prime_utrs:
                        for utr in five_prime_utrs[transcript_id]:
                            if donor_pos >= utr['start'] and donor_pos <= utr['end']:
                                splice_site_type = "5PrimeUTRSpliceSite"
                                break
                    
                    # Check if donor site is in 3' UTR if not found in 5' UTR
                    if splice_site_type == "SpliceSite" and transcript_id in three_prime_utrs:
                        for utr in three_prime_utrs[transcript_id]:
                            if donor_pos >= utr['start'] and donor_pos <= utr['end']:
                                splice_site_type = "3PrimeUTRSpliceSite"
                                break
                
                # Add donor splice site with appropriate label
                if is_lnc_RNA:
                    elements.append([chrom, donor_pos, donor_pos + 2, strand, f"{transcript_id}-lncRNAspliceSite"])
                else:
                    elements.append([chrom, donor_pos, donor_pos + 2, strand, f"{transcript_id}-{splice_site_type}"])
                
                # Acceptor site (start of next exon)
                acceptor_pos = next_exon['start']
                splice_site_type = "SpliceSite"  # Reset label
                
                if not is_lnc_RNA:  # Only check UTRs for protein-coding genes
                    # Check if acceptor site is in 5' UTR
                    if transcript_id in five_prime_utrs:
                        for utr in five_prime_utrs[transcript_id]:
                            if acceptor_pos >= utr['start'] and acceptor_pos <= utr['end']:
                                splice_site_type = "5PrimeUTRSpliceSite"
                                break
                    
                    # Check if acceptor site is in 3' UTR if not found in 5' UTR
                    if splice_site_type == "SpliceSite" and transcript_id in three_prime_utrs:
                        for utr in three_prime_utrs[transcript_id]:
                            if acceptor_pos >= utr['start'] and acceptor_pos <= utr['end']:
                                splice_site_type = "3PrimeUTRSpliceSite"
                                break
                
                # Add acceptor splice site with appropriate label
                if is_lnc_RNA:
                    elements.append([chrom, acceptor_pos - 3, acceptor_pos - 1, strand, f"{transcript_id}-lncRNAspliceSite"])
                else:
                    elements.append([chrom, acceptor_pos - 3, acceptor_pos - 1, strand, f"{transcript_id}-{splice_site_type}"])

    # Process other feature types
    with open(gff_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                fields = line.split('\t')
                if fields[2] in ['miRNA', 'scRNA', 'snRNA', 'snoRNA', 'CDS', 
                               'five_prime_UTR', 'three_prime_UTR', 'exon']:
                    chrom = validate_and_fix_contig(fields[0])
                    if chrom is None:
                        continue

                    transcript_id = next((item.split(':')[1] for item in fields[8].split(';') 
                                        if 'transcript:' in item), '').strip().split('.')[0]
                    
                    if not transcript_id:
                        continue

                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]

                    # Handle lncRNA exons - use consistent 'lncRNAexon' label
                    if fields[2] == 'exon' and transcript_id in selected_lnc_RNAs:
                        elements.append([chrom, start - 1, end, strand, f"{transcript_id}-lncRNAexon"])
                        continue

                    # Handle other features (excluding lncRNA-related features)
                    if transcript_id not in selected_lnc_RNAs:
                        region_type = ('CDS' if fields[2] == 'CDS' else 
                                    '5PrimeUTR' if fields[2] == 'five_prime_UTR' else 
                                    '3PrimeUTR' if fields[2] == 'three_prime_UTR' else fields[2])

                        if region_type in ['CDS', '5PrimeUTR', '3PrimeUTR'] and transcript_id not in mane_transcripts:
                            continue

                        if region_type != 'exon':  # Only add non-exon features for non-lncRNA transcripts
                            elements.append([chrom, start - 1, end, strand, f"{transcript_id}-{region_type}"])

    return elements

def parse_regulatory_gff(regulatory_gff_path):
    """Parse regulatory features from Ensembl GFF3 file."""
    elements = []
    
    with open(regulatory_gff_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                fields = line.split('\t')
                if fields[2] in ['promoter', 'enhancer']:
                    chrom = validate_and_fix_contig(fields[0])
                    if chrom is None:
                        continue
                    
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    
                    # Extract feature ID and extended start/end from attributes
                    attributes = {attr.split('=')[0]: attr.split('=')[1] for attr in fields[8].split(';')}
                    feature_id = attributes.get('ID', 'unknown')
                    extended_start = int(attributes.get('extended_start', start))
                    extended_end = int(attributes.get('extended_end', end))
                    
                    # Add feature annotation with appropriate type
                    feature_type = fields[2].capitalize()  # 'promoter' -> 'Promoter', 'enhancer' -> 'Enhancer'
                    elements.append([chrom, start - 1, end, strand, f"{feature_id}-{feature_type}"])
                    
                    # Adjust extended region to exclude the core region
                    if extended_start < start:
                        elements.append([chrom, extended_start - 1, start - 1, strand, 
                                      f"{feature_id}-{feature_type}Extended"])
                    if extended_end > end:
                        elements.append([chrom, end, extended_end, strand, 
                                      f"{feature_id}-{feature_type}Extended"])
    
    return elements

def count_region_types(dataframe):
    """
    Count the occurrences of each region type and distinct transcripts in the dataframe.
    """
    # Extract the region type and transcript ID from the region column
    region_types = dataframe.with_columns([
        pl.col('region').str.split('-').list.get(1).alias('region_type'),
        pl.col('region').str.split('-').list.get(0).alias('transcript_id')
    ])
    
    # Count total occurrences and distinct transcripts for each type
    counts = region_types.group_by('region_type').agg([
        pl.len().alias('total_regions'),
        pl.col('transcript_id').n_unique().alias('distinct_transcripts')
    ]).sort('total_regions', descending=True)
    
    return counts

def main(gff_path, regulatory_gff_path, output_path):
    # Parse gene annotation GFF
    regions = parse_gff(gff_path)
    
    # Parse regulatory features GFF
    if regulatory_gff_path:
        regulatory_regions = parse_regulatory_gff(regulatory_gff_path)
        regions.extend(regulatory_regions)
    
    # Create DataFrame and sort
    processed_dataframe = pl.DataFrame(
        regions, 
        schema={
            'hg38_chr': pl.Utf8, 
            'start': pl.Int64, 
            'end': pl.Int64, 
            'strand': pl.Utf8, 
            'region': pl.Utf8
        }
    ).sort(['hg38_chr', 'start', 'end'])
   
    processed_dataframe = processed_dataframe.filter(pl.col("hg38_chr") == "chr21")

    # Write output files
    processed_dataframe.write_csv(output_path, separator='\t')
    
    # Count and display region types
    region_counts = count_region_types(processed_dataframe)
    print("\n=== Region Type Counts ===")
    with pl.Config(tbl_rows=region_counts.height):
        print(region_counts)
    
    # Save counts to a file
    counts_output_path = os.path.splitext(output_path)[0] + "_counts.tsv"
    region_counts.write_csv(counts_output_path, separator='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process GFF file and create BED file')
    parser.add_argument('--gff_file', required=True, help='Path to GFF file')
    parser.add_argument('--regulatory_gff_file', help='Path to regulatory features GFF file')
    parser.add_argument('--output', required=True, help='Output BED file path')
    
    args = parser.parse_args()
    
    main(args.gff_file, args.regulatory_gff_file, args.output)
