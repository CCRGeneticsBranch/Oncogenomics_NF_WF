#!/usr/bin/env python

import pandas as pd
import glob
import sys
import os

output = sys.argv[1]
output_file = f"{output}.fusion.txt"

# Get a list of file names
file_list = glob.glob("*annotated.txt")

# Read the first file to get the header
first_file = pd.read_csv(file_list[0], delimiter="\t")
header = list(first_file.columns)

# Concatenate all files, skipping the header in all but the first file
result = pd.concat(
    [
        pd.read_csv(file, delimiter="\t", header=None, names=header, skiprows=1)
        for file in file_list
    ]
)

# Write the result to a new file
result.to_csv(output_file, index=False, sep="\t")
