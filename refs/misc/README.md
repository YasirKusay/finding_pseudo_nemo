# How to Get GNOMAD Data for IKBKG and IKBKGP1

1. Download the joint frequency chrX sites from https://gnomad.broadinstitute.org/downloads
2. Fetch both the vcf file and its index
3. Run the below command:

```
# Note, this adds some padding to IKBKG and IKBKGP1 
bcftools view -r chrX:154538238-154584046,chrX:154620391-154650275 gnomad.joint.${GNOMAD_VERSION}.sites.chrX.vcf.bgz | egrep -v "^#" > ikbkg_ikbkgp1_gnomad_entires_VCF.tsv

python3 helper_scripts/get_ikbkg_ikbkgp1_entries.py > ikbkg_ikbkgp1_gnomad_entires.tsv

rm ikbkg_ikbkgp1_gnomad_entires_VCF.tsv
```
