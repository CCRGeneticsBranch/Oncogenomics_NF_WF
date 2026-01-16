#!/usr/bin/env python3
import sys, csv, math, argparse

# Ontology terms to keep (applied ONLY to rare_annotation.tsv)
KEEP_ONTOLOGIES = {
    "frameshift_elongation",
    "frameshift_truncation",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "splice_site_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "NMD_transcript_variant",
    "stop_retained_variant",
    "lnc_RNA",
}

# Columns to remove and to use for popmax/MAF
DROP_AF_COLS = {"Global AF", "Ashkenazi Jewish AF", "Other AF"}
POP_AF_COLS = [
    "African AF",
    "East Asian AF",
    "Finnish AF",
    "Latino AF",
    "Non-Fin Eur AF",
    "South Asian AF",
]

SENTINELS = {"", ".", "NA", "N/A"}


def _to_float(x: str) -> float:
    try:
        if x is None:
            return float("nan")
        s = x.strip()
        if s in SENTINELS:
            return float("nan")
        return float(s)
    except Exception:
        return float("nan")


def extract_and_write(infile, all_outfile, rare_outfile, maf_cutoff=0.05):
    """
    Read a CRAVAT TSV, enter the `#Report level: variant` block, find the UID header,
    drop three AF columns, compute MAF across remaining population AFs,
    write full table to all_annotations, and MAF<=cutoff + ontology-filtered subset to rare_annotation.
    """
    started = False
    got_header = False
    raw_header = None
    out_header = None

    writer_all = None
    writer_rare = None

    for line in infile:
        # CRAVAT block boundary
        if line.startswith("#Report level: variant"):
            started = True
            continue
        if line.startswith("#CRAVAT Report") and started:
            break

        if started and not got_header:
            if line.startswith("UID\t"):
                got_header = True
                raw_header = [h.strip() for h in line.rstrip("\n").split("\t")]

                # Build output header: drop selected AFs, append MAF
                kept_cols = [h for h in raw_header if h not in DROP_AF_COLS]
                if "MAF" not in kept_cols:
                    kept_cols.append("MAF")
                out_header = kept_cols

                writer_all = csv.DictWriter(
                    all_outfile,
                    fieldnames=out_header,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer_rare = csv.DictWriter(
                    rare_outfile,
                    fieldnames=out_header,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer_all.writeheader()
                writer_rare.writeheader()
            continue

        # Process rows
        if got_header:
            row = line.rstrip("\n").split("\t")

            if len(row) < len(raw_header):
                row += [""] * (len(raw_header) - len(row))
            elif len(row) > len(raw_header):
                row = row[: len(raw_header)]

            rec_full = dict(zip(raw_header, row))

            # Compute MAF across remaining population AF columns
            pop_cols_present = [c for c in POP_AF_COLS if c in rec_full]
            pop_vals = [_to_float(rec_full.get(c, "")) for c in pop_cols_present]
            # Use -2 for missing MAF to indicate rare/absent from databases
            maf = max((v for v in pop_vals if not math.isnan(v)), default=-2.0)

            # Build output record: drop the three AF columns, add MAF
            out_rec = {k: rec_full.get(k, "") for k in out_header if k != "MAF"}
            out_rec["MAF"] = "" if maf == -2.0 else f"{maf:.6g}"

            # Always write full table
            writer_all.writerow(out_rec)

            # For rare table: require MAF <= cutoff AND ontology match
            # MAF of -2 indicates missing/rare variants
            ont = rec_full.get("Sequence Ontology", "").strip()
            terms = [t.strip() for t in ont.replace(";", ",").split(",") if t.strip()]
            ontology_ok = any(t in KEEP_ONTOLOGIES for t in terms)

            if ontology_ok and maf <= maf_cutoff:
                writer_rare.writerow(out_rec)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract CRAVAT variant table between report markers, "
        "drop selected AF columns, compute MAF (popmax), "
        "write full table and an ontology+MAF-filtered subset."
    )
    parser.add_argument("input", help="Input CRAVAT TSV, or '-' for stdin")
    parser.add_argument("--sample-id", help="Sample ID prefix for output files")
    parser.add_argument(
        "--all-output",
        help="Output path for full table with MAF added (default: all_annotations.txt or {sample-id}.Annotations.final.txt)",
    )
    parser.add_argument(
        "--rare-output",
        help="Output path for ontology+MAF subset (default: rare_annotation.tsv or {sample-id}_Annotations.coding.rare.txt)",
    )
    parser.add_argument(
        "--maf-cutoff",
        type=float,
        default=0.05,
        help="MAF threshold for rare output (default: 0.05)",
    )
    args = parser.parse_args()

    # Determine output file names
    if args.sample_id:
        all_default = f"{args.sample_id}.Annotations.final.txt"
        rare_default = f"{args.sample_id}.Annotations.coding.rare.txt"
    else:
        all_default = "all_annotations.txt"
        rare_default = "rare_annotation.tsv"

    all_output = args.all_output if args.all_output else all_default
    rare_output = args.rare_output if args.rare_output else rare_default

    infile = sys.stdin if args.input == "-" else open(args.input, "r")
    all_outfile = open(all_output, "w", newline="")
    rare_outfile = open(rare_output, "w", newline="")

    try:
        extract_and_write(infile, all_outfile, rare_outfile, maf_cutoff=args.maf_cutoff)
    finally:
        if infile is not sys.stdin:
            infile.close()
        all_outfile.close()
        rare_outfile.close()
