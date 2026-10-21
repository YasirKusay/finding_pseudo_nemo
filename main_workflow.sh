#!/bin/bash

set -eo pipefail
set -a
source "config.env"
set +a

# Fail early if all are missing or empty
: "${VEP_SIF:?VEP_SIF not set in config.env}"
: "${VEP_CACHE:?VEP_CACHE not set in config.env}"
: "${VEP_PLUGINS:?VEP_PLUGINS not set in config.env}"

[ $# -eq 2 ] || { echo "usage: $0 in.bam outdir" >&2; exit 1; }

infile=$1
sample=$(basename "$infile")
sample=${sample%.bam}
outfile_path=$2

if [ ! -d $outfile_path ]; then
    mkdir -p ${outfile_path}/preprocessing
else
    echo "${outfile_path} already exists. Please delete the folder or change the output directory name" >&2
    exit 1
fi

preprocessing_prefix="${outfile_path}/preprocessing/${sample}"
infile_roi_only="${preprocessing_prefix}.infile_roi_only.bam"

# also sort by name, required for safely assessing read and mate for the next step
# do not use fetch-pairs, as they could result in false positive SV calls

samtools view -h -b -M -L refs/beds/ROI.bed "$infile" -o $infile_roi_only

# samtools fastq by default ignores secondary/supplementary alignments
# while this is usually fine (primary alignment)
# if the primary part of a chimeric read is not in the roi, it will get dropped
samtools view -h $infile_roi_only | samtools sort -n | samtools fastq -1 ${preprocessing_prefix}_R1.fastq.gz \
                                                                      -2 ${preprocessing_prefix}_R2.fastq.gz \
                                                                      -s ${preprocessing_prefix}_s.fastq.gz \
                                                                      -0 ${preprocessing_prefix}_0.fastq.gz -n -

cat ${preprocessing_prefix}_s.fastq.gz ${preprocessing_prefix}_0.fastq.gz > ${preprocessing_prefix}_singletons.fastq.gz
rm ${preprocessing_prefix}_s.fastq.gz ${preprocessing_prefix}_0.fastq.gz

# realign it now to remove anything with mapq < 50

remapped_pe="${preprocessing_prefix}.remapped_PE_ONLY.bam"
remapped_singletons="${preprocessing_prefix}.remapped_singletons.bam"
remapped_merged="${preprocessing_prefix}.remapped.bam"

genome_IKBKGP1_masked="refs/indexes/genome_IKBKGP1_masked/genome_IKBKGP1_masked.fa"

# dont use -M flag as it will mark shorter split hits as secondary
# samtools depth/Vardict will drop them by default. Without it, they're
# supplementary and are kept
# -a flag will output all alternative alignments for a query sequence
# these will be flagged as secondary, so ignore it
bwa mem -Y -t 1 $genome_IKBKGP1_masked ${preprocessing_prefix}_R1.fastq.gz ${preprocessing_prefix}_R2.fastq.gz \
            | samtools view -h -b | samtools sort -@ 1 -m 4G -o $remapped_pe

# note here the flag 0 means (SE, mapped, mapped to forward strand)
# note here the flag 16 means (SE, mapped, mapped to reverse strand)
# PNEXT will always be zero
bwa mem -Y -t 1 $genome_IKBKGP1_masked ${preprocessing_prefix}_singletons.fastq.gz \
            | samtools view -h -b | samtools sort -@ 1 -m 4G -o $remapped_singletons

samtools merge -o $remapped_merged $remapped_pe $remapped_singletons
samtools index $remapped_merged

remapped_depth_file="${preprocessing_prefix}.remapped_depth.tsv"
echo -e "Chr\tPosition\tDepth" > $remapped_depth_file

# does not count duplicates by default
# does count supplementaries
samtools depth -a -b refs/beds/ROI.bed $remapped_merged >> $remapped_depth_file

# for output cleanliness
rm ${preprocessing_prefix}_R1.fastq.gz
rm ${preprocessing_prefix}_R2.fastq.gz
rm ${preprocessing_prefix}_singletons.fastq.gz
rm ${remapped_pe}
rm ${remapped_singletons}

# necessary for small variant detection, but not for CNV detection
# avoids including poorly aligned reads

outfile_prefix="${outfile_path}/${sample}"

# Note that vardict by default uses this filter flag -F 504 
# this drops unmapped, secondary and duplicate reads (if you set the -t flag, vardict will identify duplicates, that were not flagged as such)
# HOWEVER, this check will drop the supplementary reads downstream, if the -F flag is not set to 0: 
# if (!instance().conf.samfilter.equals("0") && record.getSupplementaryAlignmentFlag()) {
#     return; // Ignore the supplementary for now so that it won't skew the coverage
# }
# Hence, we manually filter it here, such that we can apply the -F flag in vardict-java and not worry about those reads

remapped_merged_filtered="${preprocessing_prefix}.remapped_filtered.bam"
samtools view -b -F 0x504 -o $remapped_merged_filtered $remapped_merged
samtools index $remapped_merged_filtered

# -c -S -E are for the bed file
# -t means dedup
# -I refers to the maximum indel length
# Combine SNVs into MNV or SNV with indel to complex if variants are located closer than `-X` option (default <=2 bases) and have good base qualities. 
# Thats why we set -X into 0, as we want to report the variants as separate records so VEP, SpliceAI, REVEL and gnomAD can annotate them individually
# currently disabled -k (local realignment) for fear of pseudogene remapping but it should be considered
# note, the bcftools norm will spit out this warning: Warning: the tag INFO/SVLEN must be defined as Number=A in order for the field to be split
# we will ignore it as vardict outputs one SV output per SV allele
# -m refers to read mismatches tolerated (greater than this, and they will be filtered out)
# we set it to 1 value higher than out cutoff in filter_poorly_aligned_reads

vardict_out="${outfile_prefix}_vardict.vcf"

# for var2vcf_valid.pl, -b command is the coordinates of the bed file it will put into the header of the vcf
# it will use the "end" as the length
vardict-java -G $genome_IKBKGP1_masked \
             -N ${sample} \
             -b $remapped_merged_filtered \
             -I 120 \
             -X 0 \
             --mfreq 0.1 \
             -t \
             --fisher \
             -f 0.10 \
             -r 5 \
             -k 0 \
             -m 9 \
             -c 1 \
             -S 2 \
             -E 3 \
             -F 0 \
             refs/beds/ikbkg_small_variant_analysis_regions.bed --nosv | \
            var2vcf_valid.pl -b refs/beds/chrX.bed -N ${sample} -A -f 0.01 | \
            bcftools sort - | bcftools norm -f $genome_IKBKGP1_masked -m - -o ${vardict_out}

# follow guide here to manually install vep caches: https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
# only use refseq at this stage
# need to also fetch Homo_sapiens.GRCh38.dna.toplevel.fa.gz from https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/
# put it in ~/vep_data/homo_sapiens_[merged|refseq]/115_GRCh38/.
# removed merged flag, only want refseq with refseq select

rm $remapped_merged_filtered
rm "${remapped_merged_filtered}.bai"

vardict_annotated_temp1="${outfile_prefix}_vardict_temp1.vcf"
vardict_annotated_temp2="${outfile_prefix}_vardict_temp2.vcf"
vardict_annotated_temp3="${outfile_prefix}_vardict_temp3.vcf"
vardict_annotated_filtered="${outfile_prefix}_vardict_annotated_filtered.vcf"
vardict_annotated_filtered_tsv="${outfile_prefix}_vardict_annotated_filtered.tsv"

singularity exec $VEP_SIF vep \
    --dir $VEP_CACHE \
    --cache \
    --offline \
    --format vcf \
    --vcf \
    --force_overwrite \
    --input_file ${vardict_out} \
    --output_file ${vardict_annotated_temp1} \
    --refseq \
    --mane \
    --species homo_sapiens \
    --hgvs \
    --symbol \
    --chr chrX \
    --assembly GRCh38 \
    --plugin NMD \
    --plugin SpliceAI,snv=$VEP_PLUGINS/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$VEP_PLUGINS/spliceai_scores.raw.indel.hg38.vcf.gz \
    --plugin AlphaMissense,file=$VEP_PLUGINS/AlphaMissense_hg38.tsv.gz,cols=all \
    --plugin REVEL,file=$VEP_PLUGINS/revel/new_tabbed_revel_grch38.tsv.gz

# MANE_Select transcript is NM_001099857.5
# Just keep the MANE_Select transcript
singularity exec $VEP_SIF filter_vep \
    -i ${vardict_annotated_temp1} \
    -o ${vardict_annotated_temp2} \
    -filter "Feature is NM_001099857.5" \
    --only_matched

bcftools +split-vep ${vardict_annotated_temp2} \
  -c Consequence,IMPACT,HGVSc,HGVSp,REVEL,am_pathogenicity,am_class,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL \
  -s mane \
  -o ${vardict_annotated_temp3}

# third flag is the gnomad threshold
python3 helper_scripts/remove_common_variants.py ${vardict_annotated_temp3} ${vardict_annotated_filtered} 0.01

# have a tsv option for the output as well
python3 convert_annotation_vcf_to_tsv.py ${vardict_annotated_filtered} ${vardict_annotated_filtered_tsv}

# clean up
rm $vardict_annotated_temp1 
rm $vardict_annotated_temp2 
rm $vardict_annotated_temp3
