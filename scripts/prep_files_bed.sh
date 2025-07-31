#!/bin/bash

aria2c https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gff3.gz
aria2c https://ftp.ensembl.org/pub/release-113/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v113.gff3.gz
gunzip Homo_sapiens.GRCh38.113.chr.gff3.gz
gunzip Homo_sapiens.GRCh38.regulatory_features.v113.gff3.gz
