import sys
from cyvcf2 import VCF, Writer

input_vcf = sys.argv[1]
output_vcf = sys.argv[2]
threshold = float(sys.argv[3])

common_variants_list = []
with open("refs/misc/common_variants_all_workflows.txt", "r") as f:
    for line in f:
        common_variants_list.append(line.strip())

with open("refs/misc/ikbkg_ikbkgp1_gnomad_entires.tsv", "r") as f:
    f.readline()
    for line in f:
        curr = line.strip().split("\t")
        curr_v = "_".join([curr[0], curr[1], curr[2], curr[3]])
        if float(curr[7]) >= threshold:
            if curr_v not in common_variants_list:
                common_variants_list.append(curr_v)

vcf = VCF(input_vcf)
writer = Writer(output_vcf, vcf)

for v in vcf:
    # assume that multi-allelic split has already happeneds
    curr_variant = (v.CHROM, str(v.POS), v.REF, v.ALT[0])
    curr_variant_s = "_".join(curr_variant)
    if curr_variant_s not in common_variants_list:
        writer.write_record(v)

writer.close()
vcf.close()