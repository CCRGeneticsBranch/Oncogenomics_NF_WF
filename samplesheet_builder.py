#!/bin/env python3
import os
import csv
from collections import defaultdict
import sys
import re

# Set default directories
DEFAULT_SAMPLESHEET_DIR = "/data/khanlab/projects/DATA/Sequencing_Tracking_Master"
DEFAULT_INPUT_DIR = "/data/khanlab/projects/DATA"
# DEFAULT_INPUT_DIR = "/data/khanlab2/kids_first_RMS/cavatica/exome_bam/DATA"

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} <patient_id> <case_name>")
    print(f"Default Mastersheet Directory: {DEFAULT_SAMPLESHEET_DIR}")
    print(f"Default Input Directory: {DEFAULT_INPUT_DIR}")
    print("To use custom directories, modify the script:")
    print(f"   - Change 'DEFAULT_SAMPLESHEET_DIR' to your samplesheet directory path")
    print(f"   - Change 'DEFAULT_INPUT_DIR' to your input directory path")
    sys.exit(1)

# Extract arguments from command-line
sample_id = sys.argv[1]
case_name = sys.argv[2]

# Use default directories
samplesheet_dir = DEFAULT_SAMPLESHEET_DIR
inputdir = DEFAULT_INPUT_DIR

# Print debug information to verify arguments and directories
print(f"Sample ID: {sample_id}")
print(f"Case Name: {case_name}")
print(f"Samplesheet Directory: {samplesheet_dir}")
print(f"Input Directory: {inputdir}")


def read_and_map_samplesheet(
    samplesheet, inputdir, column_mapping, sample_id, case_name
):
    samplesheet_data = []
    invalid_paths = []
    try:
        with open(
            samplesheet, "r", newline="", encoding="utf-8", errors="replace"
        ) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                # column mapping
                mapped_row = {
                    new_column: row.get(old_column, "").strip()
                    for old_column, new_column in column_mapping.items()
                }
                if "type" in mapped_row:
                    mapped_row["type"] = mapped_row["type"].replace(" ", "_")
                if "type" in mapped_row and mapped_row["type"] == "blood_DNA":
                    mapped_row["type"] = "normal_DNA"
                if "Diagnosis" in mapped_row:
                    mapped_row["Diagnosis"] = re.sub(
                        r",\s*", "_", mapped_row["Diagnosis"]
                    )
                    mapped_row["Diagnosis"] = mapped_row["Diagnosis"].replace(" ", ".")

                # Check if casename column has comma-separated values
                if "," in mapped_row.get("casename", ""):
                    casenames = mapped_row["casename"].split(",")
                    for individual_casename in casenames:
                        new_row = dict(mapped_row)  # Create copy of the original row
                        new_row[
                            "casename"
                        ] = individual_casename.strip()  # Update casename
                        samplesheet_data.append(new_row)
                else:
                    samplesheet_data.append(mapped_row)

                if "/" in mapped_row.get("seq_type", ""):
                    seqtypes = mapped_row["seq_type"].split("/")
                    for individual_seqtype in seqtypes[1:]:
                        new_row = dict(mapped_row)  # Create copy of the original row
                        new_row[
                            "seq_type"
                        ] = individual_seqtype.strip()  # Update seqtype
                        samplesheet_data.append(new_row)
                else:
                    samplesheet_data.append(mapped_row)

    except UnicodeDecodeError:
        print(f"Error decoding file: {samplesheet}. Please check the file encoding.")
        return [], []
    ALLOWED_SAMPLE_CAPTURES = {
        "access",
        "polya_stranded",
        "polya",
        "ribozero",
        "SmartRNA",
        "ribodepleted_nebnext_v2",
        "clin.ex.v1",
        "seqcapez.hu.ex.v3",
        "seqcapez.rms.v1",
        "agilent.v7",
        "idt_v2_plus",
        "xgen-hyb-panelv2",
        "comp_ex_v1",
        "seqcapez.hu.ex.utr.v1",
    }
    # Filter rows matching sample_id and case_name
    filtered_samplesheet_data = []
    for row in samplesheet_data:
        # if row.get("sample") == sample_id and row.get("casename") == case_name:
        if (
            row.get("sample") == sample_id
            and row.get("casename") == case_name
            and row.get("seq_type") in ["E-il", "P-il", "T-il"]
        ):
            sample_captures = row.get("sample_captures", "").strip()

            # Check if sample_captures is missing or empty
            if not sample_captures:
                print(
                    f"ERROR: 'sample_captures' is missing or empty for sample {sample_id}, case {case_name}"
                )
                sys.exit(1)
            if sample_captures not in ALLOWED_SAMPLE_CAPTURES:
                print(
                    f"ERROR: Unrecognized 'sample_captures' value '{sample_captures}' for sample {sample_id}, case {case_name}."
                )
                print(
                    "Please add this sample capture type to the workflow before creating the samplesheet."
                )
                sys.exit(1)
            if row.get("genome") not in ["hg19", "hg38", "mm10"]:
                print("genome information missing, defaulting to hg19")
                row["genome"] = "hg19"
            library_id = row["library"]
            fcid = row.get("FCID", "")
            if fcid:  # If FCID is not empty
                read1 = f"{inputdir}/Sample_{library_id}_{fcid}/Sample_{library_id}_{fcid}_R1.fastq.gz"
                read2 = f"{inputdir}/Sample_{library_id}_{fcid}/Sample_{library_id}_{fcid}_R2.fastq.gz"
            else:  # If FCID is empty
                read1 = (
                    f"{inputdir}/Sample_{library_id}/Sample_{library_id}_R1.fastq.gz"
                )
                read2 = (
                    f"{inputdir}/Sample_{library_id}/Sample_{library_id}_R2.fastq.gz"
                )

            # Check if input fastq path exists
            if os.path.exists(read1) and os.path.exists(read2):
                row["read1"] = read1
                row["read2"] = read2
                filtered_samplesheet_data.append(row)
            else:
                invalid_paths.append((read1, read2))

    return filtered_samplesheet_data, invalid_paths


def process_samplesheets(
    samplesheet_dir, inputdir, column_mapping, sample_id, case_name
):
    all_filtered_data = []
    all_invalid_paths = []

    for filename in os.listdir(samplesheet_dir):
        filepath = os.path.join(samplesheet_dir, filename)
        if os.path.isfile(filepath):
            print(f"Processing file: {filepath}")
            filtered_data, invalid_paths = read_and_map_samplesheet(
                filepath, inputdir, column_mapping, sample_id, case_name
            )
            all_filtered_data.extend(filtered_data)
            all_invalid_paths.extend(invalid_paths)

    # Group data and write to output files
    grouped_data = defaultdict(list)
    for row in all_filtered_data:
        grouped_data[(row["sample"], row["casename"])].append(row)

    for key in grouped_data:
        grouped_data[key] = list(
            {tuple(d.items()): d for d in grouped_data[key]}.values()
        )

    if not all_invalid_paths:  # Only proceed if no invalid paths were found
        for sample_casename, rows in grouped_data.items():
            if sample_casename[1].startswith("patient_"):
                if sample_casename in all_invalid_paths:
                    continue
            output_file = os.path.join(
                os.getcwd(), f"{sample_casename[0]}_{sample_casename[1]}.csv"
            )
            with open(output_file, "w", newline="") as file:
                writer = csv.DictWriter(file, fieldnames=column_mapping.values())
                writer.writeheader()
                writer.writerows(rows)
    else:
        print("Skipping CSV generation due to invalid FASTQ paths.")

    # Report invalid paths
    if all_invalid_paths:
        print("The following fastq file paths are invalid:")
        for read1, read2 in all_invalid_paths:
            print(f"Read1: {read1}, Read2: {read2}")


# column mapping
column_mapping = {
    "Patient ID": "sample",
    "Library ID": "library",
    "read1": "read1",
    "read2": "read2",
    "Enrichment step": "sample_captures",
    "Matched RNA-seq lib": "Matched_RNA",
    "Matched normal": "Matched_normal",
    "Diagnosis": "Diagnosis",
    "Case Name": "casename",
    "Type": "type",
    "FCID": "FCID",
    "Type of sequencing": "seq_type",
    "SampleRef": "genome",
}

# Call the function to process all samplesheets
process_samplesheets(samplesheet_dir, inputdir, column_mapping, sample_id, case_name)
