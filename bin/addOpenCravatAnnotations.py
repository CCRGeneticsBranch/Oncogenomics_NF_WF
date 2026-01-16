#!/usr/bin/env python3
"""
Script to merge OpenCRAVAT filtered annotations with SNPeff txt file
Replicates the logic of addAnnotations2vcf.pl but for OpenCRAVAT output

Usage: addOpenCravatAnnotations.py <opencravat_filtered.txt> <snpeff.txt>
Output: prints to stdout (redirect to file)
"""

import sys


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: addOpenCravatAnnotations.py <opencravat_filtered.txt> <snpeff.txt>",
            file=sys.stderr,
        )
        sys.exit(1)

    opencravat_file = sys.argv[1]
    snpeff_file = sys.argv[2]

    # Read OpenCRAVAT annotations into a hash
    # OpenCRAVAT columns: UID, Chrom, Position, Ref Base, Alt Base, ...
    # SNPeff columns: Chr, Start, End, Ref, Alt, ...
    # Key format for matching: Chr\tStart\tEnd\tRef\tAlt
    # For OpenCRAVAT: use columns 1, 2, 2, 3, 4 (Chrom, Position, Position, Ref, Alt)
    # Value: columns 5 to end (all other annotation columns)

    annotations = {}
    opencravat_header = None

    print(f"Reading OpenCRAVAT annotations from: {opencravat_file}", file=sys.stderr)
    with open(opencravat_file, "r") as f:
        first_line = True
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")

            if len(fields) < 5:
                continue

            # Save header for later use
            if first_line:
                opencravat_header = fields
                first_line = False
                continue

            # OpenCRAVAT format: UID(0), Chrom(1), Position(2), Ref Base(3), Alt Base(4), ...
            # Create key: Chrom, Position, Position, Ref, Alt to match SNPeff format
            # (Position appears twice because SNPeff has Start, End - for SNPs they're the same)
            key = f"{fields[1]}\t{fields[2]}\t{fields[2]}\t{fields[3]}\t{fields[4]}"

            # Value: everything from column 5 onwards (all annotation columns)
            if len(fields) > 5:
                value = "\t".join(fields[5:])
            else:
                value = ""

            if key not in annotations:
                annotations[key] = value

    print(f"Loaded {len(annotations)} annotations", file=sys.stderr)

    # Read SNPeff file and merge
    # Only output lines that have matching annotations
    matched = 0
    first_line = True

    print(f"Merging with SNPeff file: {snpeff_file}", file=sys.stderr)
    with open(snpeff_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")

            if len(fields) < 5:
                continue

            # SNPeff format: Chr(0), Start(1), End(2), Ref(3), Alt(4), ...
            # Key: first 5 columns
            key = "\t".join(fields[0:5])

            # Rest of snpeff columns (from column 5 onwards)
            if len(fields) > 5:
                vcf_info = "\t".join(fields[5:])
            else:
                vcf_info = ""

            # For header line, construct from both file headers
            if first_line:
                # Construct header: Chr, Start, End, Ref, Alt + OpenCRAVAT annotation columns + SNPeff info columns
                if opencravat_header and len(opencravat_header) > 5:
                    # Use column names from OpenCRAVAT (skip UID, Chrom, Position, Ref Base, Alt Base)
                    oc_annot_header = "\t".join(opencravat_header[5:])
                else:
                    oc_annot_header = ""
                # Use column names from SNPeff header (skip Chr, Start, End, Ref, Alt)
                print(f"{key}\t{oc_annot_header}\t{vcf_info}")
                first_line = False
                continue

            # Only print if annotation exists in OpenCRAVAT file
            if key in annotations:
                print(f"{key}\t{annotations[key]}\t{vcf_info}")
                matched += 1

    print(f"Matched and output {matched} variants", file=sys.stderr)


if __name__ == "__main__":
    main()
