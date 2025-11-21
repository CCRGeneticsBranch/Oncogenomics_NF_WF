#!/usr/bin/env python3

import sys

if len(sys.argv) != 3:
    print("Usage: python merge_padded_bed.py <input.bed> <output.bed>")
    sys.exit(1)

input_bed = sys.argv[1]
output_bed = sys.argv[2]

merged_intervals = []

with open(input_bed, "r") as f:
    for line in f:
        if line.strip() == "":
            continue
        if line.startswith("chrM"):
            continue
        chrom, start, end, name, desc = line.strip().split("\t")
        start = max(0, int(start) - 20)
        end = int(end) + 20

        # If no intervals stored yet, just add the first one
        if not merged_intervals:
            merged_intervals.append([chrom, start, end, [name], [desc]])
            continue

        last = merged_intervals[-1]
        if chrom == last[0] and start <= last[2]:  # overlap
            last[2] = max(last[2], end)
            last[3].append(name)
            last[4].append(desc)
        else:
            merged_intervals.append([chrom, start, end, [name], [desc]])

# Write output
with open(output_bed, "w") as out:
    for chrom, start, end, names, descs in merged_intervals:
        out.write(f"{chrom}\t{start}\t{end}\t{'|'.join(names)}\t{'|'.join(descs)}\n")

print(f"✅ Merged BED saved to: {output_bed}")
