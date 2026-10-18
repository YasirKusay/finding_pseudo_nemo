#!/bin/bash

set -euo pipefail

index_base="refs/indexes/genome_IKBKGP1_masked"

mkdir -p ${index_base}

wget -P ${index_base} \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

hg38_path="${index_base}/GCF_000001405.40_GRCh38.p14_genomic.fna"

gunzip ${hg38_path}.gz

tail -n +2 refs/misc/GCF_000001405_40_chrom_info.tsv | awk -F'\t' '{print $9 "\t" $12}' | while IFS=$'\t' read -r nc_name chr_name; do
    escaped=${nc_name//./\\.}
    sed -i "s/$escaped/$chr_name/g" $hg38_path
done

hg38_masked_path="${index_base}/genome_IKBKGP1_masked.fa"

bedtools maskfasta -fi $hg38_path -bed refs/beds/region_to_mask.bed -fo $hg38_masked_path
rm $hg38_path

bwa index $hg38_masked_path