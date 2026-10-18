#!/usr/bin/env python3

from math import isnan

"""
Convert a VEP-annotated VCF (with cyvcf2) into a flat, human-readable TSV.

Usage:
    python convert_annotation_vcf_to_tsv.py input.vcf(.gz) output.tsv
"""

import sys
import csv


SPLICEAI_FIELDS = [
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
]

OUTPUT_COLUMNS = [
    "Chromosome",
    "Position",
    "Reference",
    "Variant",
    "Quality",
    "Filter",
    "Transcript",
    "Coding_change",
    "Protein_change",
    "Effect",
    "Impact",
    "REVEL_score",
    "AlphaMissense_score",
    "AlphaMissense_class",
    "Highest_splice_score",
]

MISSING = (None, "", ".")


def info_str(variant, key):
    """Fetch a string INFO value, treating '.'/None/'' as missing (-> '')."""
    value = variant.INFO.get(key)
    if value in MISSING:
        return ""
    return str(value)


def clean_effect(consequence):
    """'3_prime_UTR_variant' -> '3 prime UTR variant'"""
    return consequence.replace("_", " ") if consequence else ""


def highest_splice_score(variant):
    scores = []
    for field in SPLICEAI_FIELDS:
        raw = variant.INFO.get(field)
        if raw not in MISSING:
            try:
                scores.append(float(raw))
            except ValueError:
                pass
    return max(scores) if scores else ""


def sample_genotype_string(variant, sample_index):
    """Reconstruct the 'a/b' or 'a|b' GT string for one sample."""
    genotype = variant.genotypes[sample_index]
    *alleles, phased = genotype
    sep = "|" if phased else "/"
    return sep.join("." if a == -1 else str(a) for a in alleles)


def format_number(value, decimals=None):
    if value in MISSING:
        return ""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return ""
    if value == int(value):
        return str(int(value))
    if decimals is not None:
        return f"{value:.{decimals}f}"
    return str(value)


def build_sample_columns(vcf):
    columns = []
    for sample in vcf.samples:
        columns += [
            f"{sample}_genotype",
            f"{sample}_percent_reads_variant",
            f"{sample}_read_depth",
        ]
    return columns


def get_format_value(variant, field, sample_index):
    """
    Safely pull a per-sample value out of a cyvcf2 FORMAT array, which
    can come back as either a 1D or 2D numpy array depending on the
    field's declared Number in the header.
    """


    arr = variant.format(field)
    if arr is None:
        return None
    row = arr[sample_index]

    if isnan(row[0]) or row[0] == -2147483648:
        return None

    try:
        return row[0]
    except (TypeError, IndexError):
        return row


def convert(input_path, output_path):
    from cyvcf2 import VCF

    vcf = VCF(input_path)
    header = OUTPUT_COLUMNS + build_sample_columns(vcf)

    with open(output_path, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(header)

        for variant in vcf:
            print(variant.FILTER)
            row = [
                variant.CHROM,
                variant.POS,
                variant.REF,
                ",".join(variant.ALT),
                format_number(variant.QUAL, decimals=2),
                "PASS" if variant.FILTER is None else variant.FILTER,
                info_str(variant, "MANE_SELECT"),
                info_str(variant, "HGVSc"),
                info_str(variant, "HGVSp"),
                clean_effect(info_str(variant, "Consequence")),
                info_str(variant, "IMPACT"),
                format_number(variant.INFO.get("REVEL")),
                format_number(variant.INFO.get("am_pathogenicity")),
                info_str(variant, "am_class"),
                format_number(highest_splice_score(variant), decimals=2),
            ]

            for i, sample in enumerate(vcf.samples):
                genotype = sample_genotype_string(variant, i)

                jaf = get_format_value(variant, "AF", i)
                percent_reads = (
                    format_number(jaf * 100, decimals=1) if jaf is not None else ""
                )

                # print(percent_reads)

                jdp = get_format_value(variant, "DP", i)
                depth = format_number(jdp) if jdp is not None else "."

                row += [genotype, percent_reads, depth]

            writer.writerow(row)

    vcf.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} input.vcf(.gz) output.tsv", file=sys.stderr)
        sys.exit(1)

    convert(sys.argv[1], sys.argv[2])
    print(f"Wrote {sys.argv[2]}")