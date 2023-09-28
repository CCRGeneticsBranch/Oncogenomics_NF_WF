#!/usr/bin/env python
import os
import shutil
import argparse
import shutil
import sys

def split_vcf_file(input_file_path, output_folder_path):
    with open(input_file_path, 'r') as f:
        lines = f.readlines()

    num_data_lines = len(lines)  # Total number of lines

    if num_data_lines <= 400:
        # Copy the entire file to the output folder
        output_file_path = os.path.join(output_folder_path, os.path.basename(input_file_path))
        shutil.copyfile(input_file_path, output_file_path)
        return
    else:
        header_lines = [line for line in lines if line.startswith("#")]
        data_lines = [line for line in lines if not line.startswith("#")]

        num_splits = 25
        lines_per_split = num_data_lines // num_splits  # Calculate lines per split

        for i in range(num_splits):
            start_idx = i * lines_per_split
            end_idx = (i + 1) * lines_per_split
            split_data = data_lines[start_idx:end_idx]

            split_file_name = os.path.splitext(os.path.basename(input_file_path))[0] + f"_{i + 1}.vcf"
            split_file_path = os.path.join(output_folder_path, split_file_name)

            with open(split_file_path, 'w') as split_file:
                split_file.writelines(header_lines)
                split_file.writelines(split_data)

# Usage
#example: input_vcf_file = "NCI0439_T1D_E_HTNCJBGX9.final.vcf"
input_vcf_file = sys.argv[1]
output_folder = sys.argv[2]
split_vcf_file(input_vcf_file, output_folder)







