#!/bin/bash

local_tmp=$(mktemp -d ${PWD}/tmp_XXXXXX) 
export TMPDIR=${local_tmp}

# Function for Step 1: Download and Prepare Gene Annotations
step1_gene_annotations() {
    echo "Starting Step 1: Downloading and Preparing Gene Annotations"
    curl -O https://ftp.ensembl.org/pub/release-114/gff3/homo_sapiens/Homo_sapiens.GRCh38.114.gff3.gz
    gunzip Homo_sapiens.GRCh38.114.gff3.gz
    awk -F '\t' 'BEGIN {OFS="\t"} !/^#/ { if ($3 == "gene" || $3 == "ncRNA_gene" || $3 == "pseudogene" || $3 == "mRNA" || $3 == "snoRNA" || $3 == "lnc_RNA" || $3 == "pseudogenic_transcript" || $3 == "exon" || $3 == "CDS" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR") {print $1, $4-1, $5} }' Homo_sapiens.GRCh38.114.gff3 > ensembl_genes.bed
    bedtools sort -i ensembl_genes.bed | bedtools merge -i stdin > ensembl_genes_sorted_merged.bed
    echo "Step 1 Complete: ensembl_genes_sorted_merged.bed created."
}

# Function for Step 2: Mask Regulatory Elements
step2_regulatory_elements() {
    echo "Starting Step 2: Downloading and Masking Regulatory Elements"
    EMARS_URL="https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.EMARs.v114.gff.gz"
    MOTIF_FEATURES_URL="https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.motif_features.v114.gff3.gz"
    REGULATORY_FEATURES_URL="https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v114.gff3.gz"

    curl -O $EMARS_URL
    gunzip Homo_sapiens.GRCh38.EMARs.v114.gff.gz
    awk -F '\t' 'BEGIN {OFS="\t"} !/^#/ {print $1, $4-1, $5}' Homo_sapiens.GRCh38.EMARs.v114.gff > emars.bed

    curl -O $MOTIF_FEATURES_URL
    gunzip Homo_sapiens.GRCh38.motif_features.v114.gff3.gz
    awk -F '\t' 'BEGIN {OFS="\t"} !/^#/ { if ($3 == "TF_binding_site" || $3 == "binding_site" || $3 == "motif") print $1, $4-1, $5}' Homo_sapiens.GRCh38.motif_features.v114.gff3 > motif_features.bed

    curl -O $REGULATORY_FEATURES_URL
    gunzip Homo_sapiens.GRCh38.regulatory_features.v114.gff3.gz
    awk -F '\t' 'BEGIN {OFS="\t"} !/^#/ { if ($3 == "regulatory_region" || $3 == "enhancer" || $3 == "promoter" || $3 == "promoter_flanking_region" || $3 == "CTCF_binding_site" || $3 == "open_chromatin_region") print $1, $4-1, $5}' Homo_sapiens.GRCh38.regulatory_features.v114.gff3 > regulatory_features.bed

    cat emars.bed motif_features.bed regulatory_features.bed | bedtools sort -i stdin | bedtools merge -i stdin > regulatory_elements_sorted_merged.bed
    echo "Step 2 Complete: regulatory_elements_sorted_merged.bed created."
}

# Function for Step 3: Mask Conserved Elements Using 447-Way PhyloP Scores
step3_conserved_elements_phylop() {
    echo "Starting Step 3: Downloading and Masking Conserved Elements (PhyloP)"
    PHYLO_P_BW_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP447way/hg38.phyloP447way.bw"
    curl -O $PHYLO_P_BW_URL
    bigWigToBedGraph hg38.phyloP447way.bw hg38.phyloP447way.bedGraph
    awk '$4 > 2' hg38.phyloP447way.bedGraph > conserved_regions_phylop.bed
    bedtools sort -i conserved_regions_phylop.bed | bedtools merge -i stdin > conserved_regions_phylop_sorted_merged.bed
    echo "Step 3 Complete: conserved_regions_phylop_sorted_merged.bed created."
}

# Function for Step 4: Exclude Known Sequencing Artifact Regions
step4_artifact_regions() {
    echo "Starting Step 4: Downloading and Masking Known Sequencing Artifact Regions"
    ENCODE_BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
    curl -O $ENCODE_BLACKLIST_URL
    gunzip hg38-blacklist.v2.bed.gz
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' hg38-blacklist.v2.bed > hg38-blacklist.v2.3col.bed

    FIXITFELIX_URL="https://zenodo.org/records/7535298/files/FixItFelix-main.zip?download=1"
    curl -L -o FixItFelix-main.zip $FIXITFELIX_URL
    unzip FixItFelix-main.zip
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' FixItFelix-main/collapsed.bed > FixItFelix-main/collapsed.tab.bed
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' FixItFelix-main/duplicated.bed > FixItFelix-main/duplicated.tab.bed

    cat hg38-blacklist.v2.3col.bed FixItFelix-main/collapsed.tab.bed FixItFelix-main/duplicated.tab.bed | bedtools sort -i stdin | bedtools merge -i stdin > artifact_regions_sorted_merged.bed
    echo "Step 4 Complete: artifact_regions_sorted_merged.bed created."
}

# Function for Step 5: Mask Repetitive Elements
step5_repetitive_elements() {
    echo "Starting Step 5: Downloading and Masking Repetitive Elements (RepeatMasker)"
    REPEATMASKER_OUT_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz"
    curl -O $REPEATMASKER_OUT_URL
    gunzip hg38.fa.out.gz
    awk 'BEGIN{OFS="\t"} NR > 3 {print $5, $6-1, $7}' hg38.fa.out > hg38_repeatmasker_raw.bed
    bedtools sort -i hg38_repeatmasker_raw.bed | bedtools merge -i stdin > hg38_repeatmasker_sorted_merged.bed
    echo "Step 5 Complete: hg38_repeatmasker_sorted_merged.bed created."
}

# Function for Step 6: Exclude CpG Islands
step6_cpg_islands() {
    echo "Starting Step 6: Downloading and Masking CpG Islands"
    CPG_ISLAND_URL="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz"
    curl -O $CPG_ISLAND_URL
    gunzip cpgIslandExt.txt.gz
    awk 'BEGIN{OFS="\t"} {print $2, $3, $4}' cpgIslandExt.txt > hg38_cpg_islands.bed
    bedtools sort -i hg38_cpg_islands.bed | bedtools merge -i stdin > hg38_cpg_islands_sorted_merged.bed
    echo "Step 6 Complete: hg38_cpg_islands_sorted_merged.bed created."
}

# Function for Step 7: Mask Low Mappability Regions
step7_low_mappability() {
    echo "Starting Step 7: Downloading and Masking Low Mappability Regions"
    MAPPABILITY_BW_URL="https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw"
    curl -O $MAPPABILITY_BW_URL
    bigWigToBedGraph k100.Umap.MultiTrackMappability.bw hg38_mappability_k100.bedGraph
    awk 'BEGIN{OFS="\t"} $4 < 1 {print $1, $2, $3}' hg38_mappability_k100.bedGraph > hg38_low_mappability_k100.bed
    bedtools sort -i hg38_low_mappability_k100.bed | bedtools merge -i stdin > hg38_low_mappability_k100_sorted_merged.bed
    echo "Step 7 Complete: hg38_low_mappability_k100_sorted_merged.bed created."
}

# Function to create merged_roulette_nc_scaling_sites.bed
create_roulette_sites() {
    echo "Starting creation of merged_roulette_nc_scaling_sites.bed"
    wget -r -np -nH --cut-dirs=3 -P . -A "*.tsv.gz" "http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/noncoding/"

    DOWNLOAD_DIR="./noncoding"
    TMP_DIR="./roulette_processed_tmp"
    mkdir -p "${TMP_DIR}"
    FINAL_OUTPUT_BED="merged_roulette_nc_scaling_sites.bed"
    NUM_PARALLEL_JOBS=$(nproc --all 2>/dev/null || echo 128)

    export DOWNLOAD_DIR TMP_DIR

    process_single_file_roulette() {
        local gzipped_file="$1"
        local base_filename
        base_filename=$(basename "${gzipped_file}" .tsv.gz)
        local processed_tmp_file="${TMP_DIR}/${base_filename}_processed.bed"
        gunzip -c "${gzipped_file}" \
            | awk 'NR > 1 {sub(/^chr/, "", $1); printf "%s\t%s\t%s\n", $1, $2-1, $2}' \
            | sort -k1,1V -k2,2n -k3,3n \
            | uniq \
            > "${processed_tmp_file}"
        echo "${processed_tmp_file}"
    }
    export -f process_single_file_roulette

    echo "Starting parallel processing of roulette .tsv.gz files..."
    find "${DOWNLOAD_DIR}" -type f -name "chr*_*.tsv.gz" -print0 \
        | xargs -0 -P "${NUM_PARALLEL_JOBS}" -I {} bash -c 'process_single_file_roulette "$@"' _ {} \
        > "${TMP_DIR}/list_of_processed_files.txt"

    echo "All individual roulette files processed."
    echo "Merging all processed roulette files..."

    if [ ! -s "${TMP_DIR}/list_of_processed_files.txt" ]; then
        echo "Error: No processed roulette files found in ${TMP_DIR}. Exiting."
        echo "Warning: No roulette files found to process. Creating an empty ${FINAL_OUTPUT_BED}."
        touch "${FINAL_OUTPUT_BED}"
    else
        cat $(cat "${TMP_DIR}/list_of_processed_files.txt") \
            | sort -k1,1V -k2,2n -k3,3n \
            | uniq \
            | bedtools merge -i stdin \
            > "${FINAL_OUTPUT_BED}"
    fi
    echo "Creation of merged_roulette_nc_scaling_sites.bed complete."
}


# Main script execution
echo "Starting script execution. Running independent steps in parallel."

# Run steps 1-7 and roulette site creation in parallel
step1_gene_annotations &
step2_regulatory_elements &
step3_conserved_elements_phylop &
step4_artifact_regions &
step5_repetitive_elements &
step6_cpg_islands &
step7_low_mappability &
create_roulette_sites &

Wait for all background jobs to complete
echo "Waiting for all parallel data processing steps to complete..."
wait
sleep 10
echo "All parallel data processing steps finished."

# Step 8: Combine All Masks and Generate Final Neutral Regions
echo "Step 8: Combining All Masks and Generating Final Neutral Regions"
# ... existing code ...
# Download hg38 chromosome sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# Concatenate all individual mask BED files
cat \
    ensembl_genes_sorted_merged.bed \
    regulatory_elements_sorted_merged.bed \
    conserved_regions_phylop_sorted_merged.bed \
    artifact_regions_sorted_merged.bed \
    hg38_repeatmasker_sorted_merged.bed \
    hg38_cpg_islands_sorted_merged.bed \
    hg38_low_mappability_k100_sorted_merged.bed \
    > all_masks_raw.bed

awk '{sub(/^chr/, "", $1); print $0}' hg38.chrom.sizes > hg38.chrom.sizes.nochr

# Filter the combined mask to include only chromosomes present in hg38.chrom.sizes
awk '{print $1}' merged_roulette_nc_scaling_sites.bed | sort -V | uniq > valid_chroms.txt
awk 'BEGIN{OFS="\t"; while((getline < "valid_chroms.txt") > 0) valid_chrom_map[$0]=1} {current_chrom=$1; if (current_chrom in valid_chrom_map) { $1=current_chrom; print $0} }' all_masks_raw.bed > all_masks_raw_filtered.bed
awk 'BEGIN{OFS="\t"; while((getline < "valid_chroms.txt") > 0) valid_chrom_map[$0]=1} {current_chrom=$1; if (current_chrom in valid_chrom_map) { $1=current_chrom; print $0} }' hg38.chrom.sizes.nochr > hg38.chrom.sizes.nochr.sorted

# Sort (using genome file order) and merge the filtered combined mask file
# This was duplicated, keeping one.
bedtools sort -i all_masks_raw_filtered.bed -g hg38.chrom.sizes.nochr.sorted > comprehensive_mask_sorted.bed

bedtools merge -i comprehensive_mask_sorted.bed > comprehensive_mask_merged.bed

# Generate the complement of the mask to get neutral regions
bedtools complement -i comprehensive_mask_merged.bed -g hg38.chrom.sizes.nochr.sorted > neutral_regions_hg38_tmp.bed

# The following line was present earlier, it seems the roulette sites are now created by create_roulette_sites function
# bedtools sort -i roulette_nc_scaling_sites.bed | bedtools merge -i stdin > merged_roulette_nc_scaling_sites.bed
# Ensure merged_roulette_nc_scaling_sites.bed exists, created by create_roulette_sites function now.

# The duplicated sort and merge for all_masks_raw_chr_filtered.bed and complement have been removed.

## The SUBCOMMAND to create merged_roulette_nc_scaling_sites.bed is now encapsulated in the create_roulette_sites function
## and run in parallel.

# Intersect with roulette sites
# Ensure merged_roulette_nc_scaling_sites.bed is ready (it should be after 'wait')
if [ ! -f merged_roulette_nc_scaling_sites.bed ]; then
    echo "Error: merged_roulette_nc_scaling_sites.bed not found. This file should have been created by the roulette processing step."
    exit 1
fi
bedtools intersect -a neutral_regions_hg38_tmp.bed -b merged_roulette_nc_scaling_sites.bed > hq_neutral_roulette_nc_scaling_sites.bed

awk '{for(i=$2;i<$3;i++) print $1"\t"i"\t"i+1}' hq_neutral_roulette_nc_scaling_sites.bed > hq_neutral_roulette_nc_scaling_sites_exploded.bed
shuf hq_neutral_roulette_nc_scaling_sites_exploded.bed | head -n 100000000 > hq_neutral_roulette_nc_scaling_sites_exploded_sampled.bed
bedtools sort -i hq_neutral_roulette_nc_scaling_sites_exploded_sampled.bed -g hg38.chrom.sizes.nochr.sorted > hq_neutral_roulette_nc_scaling_sites_exploded_sampled_sorted.bed
bedtools merge -i hq_neutral_roulette_nc_scaling_sites_exploded_sampled_sorted.bed > hq_neutral_roulette_nc_scaling_sites_exploded_sampled_sorted_merged.bed
echo "Step 9 Complete: Add strand and region_type columns" # Renaming to Step 9 for clarity
awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3, ".", "Unannotated_NC-Unannotated_NC"}' hq_neutral_roulette_nc_scaling_sites_exploded_sampled_sorted_merged.bed > features_by_loci.bed

find . -mindepth 1 -maxdepth 1 ! -name 'features_by_loci.bed' -exec rm -rf {} +

echo "Step 8 (now Step 9 logic) Complete: Final neutral regions file features_by_loci.bed for neutral sites created." # Adjusted echo
echo "Pipeline Finished!"
