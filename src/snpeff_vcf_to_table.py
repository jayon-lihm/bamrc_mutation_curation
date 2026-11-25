#!/usr/bin/env python3
"""Flatten a SnpEff-annotated bamrc VCF into a tab-separated table.

Usage:
    snpeff_vcf_to_table.py -i Pt1.bamrc_ann.vcf -o Pt1.bamrc_ann.protein_coding.tsv [--patient Pt1]

One row per sample x alt allele x SnpEff annotation. "structural_interaction_variant",
"sequence_feature" and "intragenic_variant" annotations are skipped; an alt allele left with no other annotation still gets one
row, with empty annotation columns.

Columns:
    patient, sample, chrom, pos, ref, alt,
    mutation_type   SnpEff annotation (e.g. missense_variant, stop_gained)
    gene            SnpEff gene name
    protein_change  HGVS.p (e.g. p.G12D)
    raw_depth       total depth recorded by bam-readcount (FORMAT DP)
    depth           ref_count + alt_count
    ref_count, alt_count   from FORMAT AD for this alt allele
    vaf             alt_count / depth, 4 decimals (empty if depth is 0)

Expects the FORMAT written by convert_bamrc_to_vcf.py: ALT:GT:DP:AD[:VAF], where
AD holds one "ref,alt" pair per alt allele, separated by ";".
"""

import argparse
import csv
import os
import sys

SKIP_ANNOTATIONS = {"structural_interaction_variant", "sequence_feature", "intragenic_variant"}
COLUMNS = [
    "patient", "sample", "chrom", "pos", "ref", "alt", "mutation_type", "gene",
    "protein_change", "raw_depth", "depth", "ref_count", "alt_count", "vaf",
]


def parse_ann(info):
    """Return {alt allele: [(mutation_type, gene, protein_change), ...]} from INFO ANN."""
    ann = {}
    for field in info.split(";"):
        if not field.startswith("ANN="):
            continue
        for entry in field[4:].split(","):
            parts = entry.split("|")
            if len(parts) < 11:
                continue
            allele, effect, gene, hgvs_p = parts[0], parts[1], parts[3], parts[10]
            if set(effect.split("&")) <= SKIP_ANNOTATIONS:
                continue
            ann.setdefault(allele, []).append((effect, gene, hgvs_p))
    return ann


def sample_counts(fmt_keys, value, n_alts):
    """Return raw_depth and [(ref_count, alt_count), ...] per alt allele for one sample."""
    fields = dict(zip(fmt_keys, value.split(":")))
    dp = fields.get("DP", "0")
    raw_depth = int(dp) if dp.isdigit() else 0
    pairs = []
    for pair in fields.get("AD", "").split(";"):
        try:
            r, a = pair.split(",")
            pairs.append((int(r), int(a)))
        except ValueError:
            pairs.append((0, 0))
    pairs += [(0, 0)] * (n_alts - len(pairs))
    return raw_depth, pairs[:n_alts]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="SnpEff-annotated VCF")
    ap.add_argument("-o", "--output", required=True, help="output TSV")
    ap.add_argument("--patient", help="patient ID (default: VCF file name prefix)")
    args = ap.parse_args()

    patient = args.patient or os.path.basename(args.input).split(".")[0]
    samples = None
    n_rows = 0

    with open(args.input) as vcf, open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for line in vcf:
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = cols[9:]
                continue
            if samples is None:
                sys.exit("ERROR: no #CHROM header line before records")

            chrom, pos, _id, ref, alt_field, _q, _f, info, fmt = cols[:9]
            alts = alt_field.split(",")
            ann = parse_ann(info)
            fmt_keys = fmt.split(":")

            for sample, value in zip(samples, cols[9:]):
                raw_depth, pairs = sample_counts(fmt_keys, value, len(alts))
                for alt, (ref_count, alt_count) in zip(alts, pairs):
                    depth = ref_count + alt_count
                    vaf = f"{alt_count / depth:.4f}" if depth else ""
                    for effect, gene, hgvs_p in ann.get(alt) or [("", "", "")]:
                        writer.writerow([patient, sample, chrom, pos, ref, alt, effect, gene,
                                         hgvs_p, raw_depth, depth, ref_count, alt_count, vaf])
                        n_rows += 1

    print(f"Wrote {n_rows} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
