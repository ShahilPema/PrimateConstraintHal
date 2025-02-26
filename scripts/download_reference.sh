gca_id=$1

# Function to download a genome given a GCA ID
download_genome() {
    local gca_id=$1
    local max_retries=5
    local delay=$((RANDOM % 9 + 2))

    for ((i=1; i<=max_retries; i++)); do
        # Get the FTP path for the genome
        ftp_path=$(esearch -db assembly -query "$gca_id" | \
                   efetch -format docsum | \
                   xtract -pattern DocumentSummary -element FtpPath_GenBank | \
                   awk '{print $1; exit}')
        
        if [ -n "$ftp_path" ]; then
            break
        fi

        echo "Attempt $i failed due to rate limit or network issue. Retrying in $delay seconds..."
        sleep $delay

        # Exponential backoff
        delay=$((delay + RANDOM % 9 + 1))

        echo "Attempt $i failed. Retrying..."
        sleep 2
    done

    if [ -z "$ftp_path" ]; then
        echo "No FTP path found for $gca_id after $max_retries attempts"
        return
    fi

    # Extract the filename
    filename=$(basename "$ftp_path")_genomic.fna.gz

    # Download the genome
    wget "${ftp_path}/${filename}"

    echo "Downloaded $filename"
}

# Set NCBI API key if you have one
# export NCBI_API_KEY=your_api_key_here

echo "Processing GCA ID: $gca_id"
if [ "${gca_id}" = "GCA_007565055.1" ]; then
        #There were some contigs missing in version downloaded by this script.
        wget https://hgdownload.soe.ucsc.edu/hubs/GCF/007/565/055/GCF_007565055.1/GCF_007565055.1.fa.gz
        mv GCF_007565055.1.fa.gz GCF_007565055.1_genomic.fna.gz
elif [ "${gca_id}" = "GCA_000164805.1" ]; then
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/164/805/GCF_000164805.1_Tarsius_syrichta-2.0.1/GCF_000164805.1_Tarsius_syrichta-2.0.1_genomic.fna.gz
else
    download_genome "$gca_id"
fi
