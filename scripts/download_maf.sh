#!/bin/bash

chr=$1
url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/maf/"
file="chr${chr}.maf.gz"

# Using aria2c with optimized settings for large files
aria2c --file-allocation=none -x 16 -s 16 -j 16 -k 1M "${url}${file}"

