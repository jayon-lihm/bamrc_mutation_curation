#!/usr/bin/env python3
"""Convert one patient's bam-readcount outputs into a multi-sample VCF.

Usage:
    convert_bamrc_to_vcf.py --patient Pt1 \
        --bamrc Pt1_Tumor.bamrc.tsv Pt1_Normal.bamrc.tsv \
        [--samples Pt1_Tumor Pt1_Normal] \
        [--reference Homo_sapiens_assembly38.fasta] \
        -o Pt1.bamrc.vcf

One sample column per bamrc file, in the order given. Sample names default to
the file name without ".bamrc.tsv".

Alt alleles: every non-reference allele (excluding "=" and "N") with count > 0
in at least one of the patient's samples. Sites where all samples are 100%
reference are skipped.

FORMAT = ALT:GT:DP:AD:VAF, one entry per alt allele joined with ";":
    ALT  alt allele(s), comma-separated (GT and AD stay ";"-separated per allele)
    GT   1/0 if this sample has reads for the allele, else 0/0 (./. if the
         sample has no bamrc line at this site)
    DP   total depth (bamrc column 4)
    AD   ref_count,alt_count for that allele
    VAF  alt_count / (ref_count + alt_count) from AD, 4 decimals ("." if 0/0)

Indels (bam-readcount -i):
    "+XYZ" at POS p -> insertion after p:  REF=r, ALT=r+XYZ at POS p
    "-XYZ" at POS p -> deletion of p..:    REF=a+XYZ, ALT=a at POS p-1,
                                           where a is the reference base at p-1
                                           (needs --reference with a .fai index)
"""

import argparse
import os
import sys
from collections import OrderedDict

SKIP_ALLELES = {"=", "N"}
FORMAT = "ALT:GT:DP:AD:VAF"


class Faidx:
    """Minimal indexed FASTA reader (uses <fasta>.fai)."""

    def __init__(self, fasta):
        self.fasta = fasta
        self.index = {}
        with open(fasta + ".fai") as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self.fh = open(fasta, "rb")

    def base(self, chrom, pos):
        """1-based single reference base."""
        length, offset, linebases, linewidth = self.index[chrom]
        i = pos - 1
        self.fh.seek(offset + (i // linebases) * linewidth + i % linebases)
        return self.fh.read(1).decode().upper()


def parse_bamrc(path):
    """Return {(chrom, pos): {"ref", "depth", "counts": {allele: count}}} and chrom order."""
    sites, chroms = {}, []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, ref, depth = fields[0], int(fields[1]), fields[2].upper(), int(fields[3])
            counts = {}
            for f in fields[4:]:
                allele, count = f.split(":", 2)[:2]
                counts[allele.upper()] = int(count)
            sites[(chrom, pos)] = {"ref": ref, "depth": depth, "counts": counts}
            if chrom not in chroms:
                chroms.append(chrom)
    return sites, chroms


def sample_entry(site, ref, bamrc_alts, vcf_alts):
    """FORMAT value for one sample at one record."""
    n = len(bamrc_alts)
    if site is None:  # no bamrc line for this sample here
        gts, dp, ads, vafs = ["./."] * n, "0", ["0,0"] * n, ["."] * n
    else:
        ref_count = site["counts"].get(ref, 0)
        alt_counts = [site["counts"].get(a, 0) for a in bamrc_alts]
        gts = ["1/0" if c > 0 else "0/0" for c in alt_counts]
        dp = str(site["depth"])
        ads = [f"{ref_count},{c}" for c in alt_counts]
        vafs = [f"{c / (ref_count + c):.4f}" if ref_count + c > 0 else "." for c in alt_counts]
    return ":".join([",".join(vcf_alts), ";".join(gts), dp, ";".join(ads), ";".join(vafs)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--bamrc", nargs="+", required=True, help="bamrc files, one per sample")
    ap.add_argument("--samples", nargs="+", help="sample names (default: from file names)")
    ap.add_argument("--reference", help="reference FASTA with .fai (needed only for deletions)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    samples = args.samples or [os.path.basename(p).replace(".bamrc.tsv", "") for p in args.bamrc]
    if len(samples) != len(args.bamrc):
        sys.exit("ERROR: --samples and --bamrc must have the same length")

    per_sample, chrom_order = [], []
    for path in args.bamrc:
        sites, chroms = parse_bamrc(path)
        per_sample.append(sites)
        chrom_order += [c for c in chroms if c not in chrom_order]

    keys = sorted({k for s in per_sample for k in s}, key=lambda k: (chrom_order.index(k[0]), k[1]))
    fasta = Faidx(args.reference) if args.reference else None

    records = []
    for chrom, pos in keys:
        site_list = [s.get((chrom, pos)) for s in per_sample]
        ref = next(s["ref"] for s in site_list if s)

        # union of alt alleles with count > 0 in any sample (file order kept)
        alts = OrderedDict()
        for s in site_list:
            if not s:
                continue
            for allele, count in s["counts"].items():
                if count > 0 and allele != ref and allele not in SKIP_ALLELES:
                    alts[allele] = True

        snv_ins = [a for a in alts if not a.startswith("-")]
        dels = [a for a in alts if a.startswith("-")]

        # SNVs + insertions share POS p and REF r
        if snv_ins:
            vcf_alts = [ref + a[1:] if a.startswith("+") else a for a in snv_ins]
            records.append((chrom, pos, ref, snv_ins, vcf_alts, site_list))

        # deletions anchored at p-1
        if dels:
            if fasta is None:
                sys.exit(f"ERROR: deletion at {chrom}:{pos} needs --reference")
            anchor = fasta.base(chrom, pos - 1)
            longest = max((d[1:] for d in dels), key=len)
            del_ref = anchor + longest
            vcf_alts = [anchor + longest[len(d) - 1:] for d in dels]
            records.append((chrom, pos - 1, del_ref, dels, vcf_alts, site_list, ref))

    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=convert_bamrc_to_vcf.py patient={args.patient}\n")
        out.write('##FORMAT=<ID=ALT,Number=.,Type=String,Description="Alt allele(s), comma-separated">\n')
        out.write('##FORMAT=<ID=GT,Number=.,Type=String,Description="1/0 if this sample has reads for the alt allele, else 0/0; \';\'-separated per allele">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth (bam-readcount column 4)">\n')
        out.write('##FORMAT=<ID=AD,Number=.,Type=String,Description="ref_count,alt_count per alt allele, \';\'-separated">\n')
        out.write('##FORMAT=<ID=VAF,Number=.,Type=String,Description="alt/(ref+alt) from AD per alt allele, \';\'-separated">\n')
        out.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples) + "\n")
        records.sort(key=lambda r: (chrom_order.index(r[0]), r[1]))
        for rec in records:
            chrom, vpos, vref, bamrc_alts, vcf_alts, site_list = rec[:6]
            count_ref = rec[6] if len(rec) > 6 else vref  # ref base in bamrc counts
            vid = ";".join(f"{chrom}_{vpos}_{vref}_{a}" for a in vcf_alts)
            cols = [chrom, str(vpos), vid, vref, ",".join(vcf_alts), ".", ".", ".", FORMAT]
            cols += [sample_entry(s, count_ref, bamrc_alts, vcf_alts) for s in site_list]
            out.write("\t".join(cols) + "\n")


if __name__ == "__main__":
    main()
