#!/usr/bin/python

import os
import pysam
import re
import shlex
import string
import subprocess
import tempfile
import pandas as pd
import sys
from numpy import array
from sys import argv


hsmetrics = sys.argv[7]
df = pd.read_csv(hsmetrics, comment="#", sep="\t")

# Extract column values for MEAN_BAIT_COVERAGE and MEAN_TARGET_COVERAGE
mean_bait_coverage = df["MEAN_BAIT_COVERAGE"].values
mean_target_coverage = df["MEAN_TARGET_COVERAGE"].values

# Round the values
rounded_mean_bait_coverage = mean_bait_coverage.round().astype(int)
rounded_mean_target_coverage = mean_target_coverage.round().astype(int)

# Convert rounded values to strings
rounded_mean_bait_str = rounded_mean_bait_coverage.astype(str)
rounded_mean_target_str = rounded_mean_target_coverage.astype(str)


def _count_reads(bam):
    stats = {"total_reads": 0, "mapped_reads": 0}

    with pysam.Samfile(bam, "rb") as bam_file:
        for read in bam_file:
            stats["total_reads"] += 1
            if not read.is_unmapped:
                stats["mapped_reads"] += 1

    stats["percent_mapped"] = (
        float(stats["mapped_reads"]) / float(stats["total_reads"])
    ) * 100
    stats["percent_mapped"] = "{:.2f}".format(stats["percent_mapped"])

    return stats


def _coverage(bam):
    ##samtools     = '/usr/local/apps/samtools/1.2/bin/samtools';
    samtools = "samtools"
    ##intersectbed = '/usr/local/apps/bedtools/2.24.0/bin/intersectBed'
    intersectbed = "bedtools intersect"
    targets = argv[2]

    cmd = f'{samtools} view -F 0x0400 -b "{bam}" | {samtools} depth -q 20 -Q 30 -b {targets} /dev/stdin | awk -v OFS="\\t" \'{{print $1,$2-1,$2,$3}}\' | {intersectbed} -a stdin -b {targets} -wo'

    # cmd = '%s view -F 0x0400 -b "%s" | %s depth -q 20 -Q 30 -b %s /dev/stdin | awk -v OFS="\\t" \'{print $1,$2-1,$2,$3}\' | %s -a stdin -b %s -wo'  % (samtools, bam, samtools, targets, intersectbed, targets)

    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    # proc = subprocess.Popen(["-c", cmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    stats = {}

    total_positions = 0
    positions_at_1000 = 0
    positions_at_400 = 0
    positions_at_200 = 0
    positions_at_100 = 0
    positions_at_50 = 0
    positions_at_30 = 0
    positions_at_20 = 0
    positions_at_15 = 0
    positions_at_10 = 0
    positions_at_5 = 0

    targets_fh = open(targets, "r")

    for line in targets_fh:
        line = line.rstrip(os.linesep)

        cols = line.split("\t")

        start = int(cols[1])
        end = int(cols[2])

        total_positions += end - start

    for line in proc.stdout:
        line = line.decode().rstrip(os.linesep)

        cols = line.split("\t")

        chr = cols[0]
        pos = cols[2]
        depth = int(cols[3])

        if depth >= 1000:
            positions_at_1000 += 1

        if depth >= 400:
            positions_at_400 += 1

        if depth >= 200:
            positions_at_200 += 1

        if depth >= 100:
            positions_at_100 += 1

        if depth >= 50:
            positions_at_50 += 1

        if depth >= 30:
            positions_at_30 += 1

        if depth >= 20:
            positions_at_20 += 1

        if depth >= 15:
            positions_at_15 += 1

        if depth >= 10:
            positions_at_10 += 1

        if depth >= 5:
            positions_at_5 += 1

    stats["percent_at_1000"] = (float(positions_at_1000) / float(total_positions)) * 100
    stats["percent_at_1000"] = "{:.2f}".format(stats["percent_at_1000"])

    stats["percent_at_400"] = (float(positions_at_400) / float(total_positions)) * 100
    stats["percent_at_400"] = "{:.2f}".format(stats["percent_at_400"])

    stats["percent_at_200"] = (float(positions_at_200) / float(total_positions)) * 100
    stats["percent_at_200"] = "{:.2f}".format(stats["percent_at_200"])

    stats["percent_at_100"] = (float(positions_at_100) / float(total_positions)) * 100
    stats["percent_at_100"] = "{:.2f}".format(stats["percent_at_100"])

    stats["percent_at_50"] = (float(positions_at_50) / float(total_positions)) * 100
    stats["percent_at_50"] = "{:.2f}".format(stats["percent_at_50"])

    stats["percent_at_30"] = (float(positions_at_30) / float(total_positions)) * 100
    stats["percent_at_30"] = "{:.2f}".format(stats["percent_at_30"])

    stats["percent_at_20"] = (float(positions_at_20) / float(total_positions)) * 100
    stats["percent_at_20"] = "{:.2f}".format(stats["percent_at_20"])

    stats["percent_at_15"] = (float(positions_at_15) / float(total_positions)) * 100
    stats["percent_at_15"] = "{:.2f}".format(stats["percent_at_15"])

    stats["percent_at_10"] = (float(positions_at_10) / float(total_positions)) * 100
    stats["percent_at_10"] = "{:.2f}".format(stats["percent_at_10"])

    stats["percent_at_5"] = (float(positions_at_5) / float(total_positions)) * 100
    stats["percent_at_5"] = "{:.2f}".format(stats["percent_at_5"])

    return stats


def _ontarget(bam):
    ##intersectbed = '/usr/local/apps/bedtools/2.24.0/bin/intersectBed'
    intersectbed = "bedtools intersect"
    targets = argv[2]
    ontarget_bam = "%s_ontarget.bam" % (os.path.basename(bam))

    ## os.getenv('USER', default_value)
    temp = tempfile.NamedTemporaryFile(dir=argv[3])
    cmd = '%s -abam "%s" -b %s' % (intersectbed, bam, targets)
    args = shlex.split(cmd)

    subprocess.check_call(args, stdout=temp)

    temp.seek(0)

    stats = {
        "ontarget_reads": 0,
        "unique_ontarget_reads": 0,
        "hq_unique_ontarget_reads": 0,
    }

    min_mapq = 999999999
    mapq_array = []

    with pysam.Samfile(temp.name, "rb") as bam_file:
        for read in bam_file:
            stats["ontarget_reads"] += 1
            if not read.is_duplicate:
                stats["unique_ontarget_reads"] += 1
                mapq_array.append(read.mapq)
                if read.mapq < min_mapq:
                    min_mapq = read.mapq
                if read.mapq >= 30:
                    stats["hq_unique_ontarget_reads"] += 1

    numpy_array = array(mapq_array)

    stats["min_mapq"] = min_mapq
    stats["mean_mapq"] = numpy_array.mean()
    stats["mean_mapq"] = "{:.2f}".format(stats["mean_mapq"])
    stats["stddev_mapq"] = numpy_array.std()
    stats["stddev_mapq"] = "{:.2f}".format(stats["stddev_mapq"])

    return stats


print(
    "\t".join(
        [
            "#Patient",
            "Sample",
            "Diagnosis",
            "MEAN_BAIT_COVERAGE",
            "MEAN_TARGET_COVERAGE",
            "total_reads",
            "mapped_reads",
            "percent_mapped",
            "ontarget_reads",
            "percent_ontarget",
            "unique_ontarget_reads",
            "percent_unique_ontarget",
            "min_mapq",
            "mean_mapq",
            "stddev_mapq",
            "hq_unique_ontarget_reads",
            "percent_hq_unique_ontarget",
            "percent_hq_unique_positions_at_5x",
            "percent_hq_unique_positions_at_10x",
            "percent_hq_unique_positions_at_15x",
            "percent_hq_unique_positions_at_20x",
            "percent_hq_unique_positions_at_30x",
            "percent_hq_unique_positions_at_50x",
            "percent_hq_unique_positions_at_100x",
            "percent_hq_unique_positions_at_200x",
            "percent_hq_unique_positions_at_400x",
            "percent_hq_unique_positions_at_1000x",
        ]
    )
)

for filename in [argv[1]]:
    m = re.search("(.+)\.final\.bam$", filename)

    if m is not None:
        sample = m.group(1)
        stats1 = _count_reads(filename)
        stats2 = _ontarget(filename)
        stats3 = _coverage(filename)

        stats2["percent_ontarget"] = (
            float(stats2["ontarget_reads"]) / float(stats1["mapped_reads"])
        ) * 100
        stats2["percent_ontarget"] = "{:.2f}".format(stats2["percent_ontarget"])

        stats2["percent_unique_ontarget"] = (
            float(stats2["unique_ontarget_reads"]) / float(stats2["ontarget_reads"])
        ) * 100
        stats2["percent_unique_ontarget"] = "{:.2f}".format(
            stats2["percent_unique_ontarget"]
        )

        stats2["percent_hq_unique_ontarget"] = (
            float(stats2["hq_unique_ontarget_reads"])
            / float(stats2["unique_ontarget_reads"])
        ) * 100
        stats2["percent_hq_unique_ontarget"] = "{:.2f}".format(
            stats2["percent_hq_unique_ontarget"]
        )

        print(
            "\t".join(
                [
                    argv[4],
                    argv[5],
                    argv[6],
                    rounded_mean_bait_str[0],
                    rounded_mean_target_str[0],
                    str(stats1["total_reads"]),
                    str(stats1["mapped_reads"]),
                    str(stats1["percent_mapped"]),
                    str(stats2["ontarget_reads"]),
                    str(stats2["percent_ontarget"]),
                    str(stats2["unique_ontarget_reads"]),
                    str(stats2["percent_unique_ontarget"]),
                    str(stats2["min_mapq"]),
                    str(stats2["mean_mapq"]),
                    str(stats2["stddev_mapq"]),
                    str(stats2["hq_unique_ontarget_reads"]),
                    str(stats2["percent_hq_unique_ontarget"]),
                    str(stats3["percent_at_5"]),
                    str(stats3["percent_at_10"]),
                    str(stats3["percent_at_15"]),
                    str(stats3["percent_at_20"]),
                    str(stats3["percent_at_30"]),
                    str(stats3["percent_at_50"]),
                    str(stats3["percent_at_100"]),
                    str(stats3["percent_at_200"]),
                    str(stats3["percent_at_400"]),
                    str(stats3["percent_at_1000"]),
                ]
            )
        )
