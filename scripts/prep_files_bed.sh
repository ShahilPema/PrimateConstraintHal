#!/bin/bash

aria2c https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.summary.txt.gz
zcat < MANE.GRCh38.v1.3.summary.txt.gz | cut -f8 | tail -n +2 | cut -d"." -f1 > MANE_transcript.txt
aria2c https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.chr.gff3.gz
gunzip Homo_sapiens.GRCh38.112.chr.gff3.gz
