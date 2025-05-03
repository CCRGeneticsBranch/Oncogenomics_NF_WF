#!/usr/bin/env python
"""
Script imported from COMPASS team pipeline https://github.com/manGitHub/ngs_pipeline_COMPASS_master/blob/master/scripts/cnvkit_geneann_ngs.py
> $Pipeline/CP12361/MD-23-6345/CP12361_T2D_E/sequenza/CP12361_T2D_E/CP12361_T2D_E_segments.txt \
> /data/Compass/Ref/hg19/design/NextraExome_v1.target.hg19.merged.bed \
> /data/Compass/Ref/hg19/annotation/geneLists/combinedList_05312018 \
> output/CP12361_T2D_E_temp1.txt \
> output/CP12361_T2D_E_temp2.txt \
> output/CP12361_T2D_E_temp3.txt \
> output/CP12361_T2D_E_temp4.txt \
> output/CP12361_T2D_E_temp5.txt \
> output/CP12361_T2D_E_temp6.txt \
> output/CP12361_T2D_E.txt \
> output/CP12361_T2D_E_genelevel.txt
"""

import sys
from sys import argv
import os
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
import re
import contextlib
import io
import warnings

warnings.filterwarnings("ignore", module="numpy")


# Reading all the input
lib_id = sys.argv[1]
output_dir = sys.argv[2]
# input_file = sys.argv[1] + '.bed'
bed_file_input = sys.argv[3]
gene_list = sys.argv[4]

""" New Style Call
#Reading all the input files:
cns = "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E.call.cns"
#input_file = "CP11323_T1D_E.call.bed"
bed_file_input= "/data/Compass/Ref/hg19/design/NextraExome_v1.target.hg19.merged.bed"
gene_list = "/data/Compass/Ref/hg19/annotation/geneLists/combinedList_05312018"
output_dir = "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/"
lib_id = "CP11323_T1D_E"

"""

# assigning all the temp files:
# cns = f'{output_dir}/{lib_id}.call.cns'
cns = sys.argv[5]
segment_temp1 = f"{output_dir}/{lib_id}_temp1.txt"
py_intersection1 = f"{output_dir}/{lib_id}_temp2.txt"
result_temp1 = f"{output_dir}/{lib_id}_temp3.txt"
gene_seperated_file = f"{output_dir}/{lib_id}_temp4.txt"
py_intersection2 = f"{output_dir}/{lib_id}_temp5.txt"
final_temp = f"{output_dir}/{lib_id}_temp6.txt"
output_txt = f"{output_dir}/{lib_id}.txt"
output_gl = f"{output_dir}/{lib_id}_genelevel.txt"
cns_temp1 = f"{output_dir}/cns.tmp"

""" Old style call
#Reading all the temp files:
segment_temp1= "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_temp1.txt"
py_intersection1= "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_temp2.txt"
result_temp1= "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_temp3.txt"
gene_seperated_file= "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_temp4.txt"
py_intersection2 = "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_temp5.txt"
final_temp= "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_temp6.txt"
output_txt = "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E.txt"
output_gl = "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/CP11323_T1D_E_genelevel.txt"
#cns = "CP11323_T1D_E.call.cns"
cns_temp1 = "/data/Compass/dev/valdezkm/ngs/CP11323/MD-24-2066/CP11323_T1D_E/cnvkit/cns.tmp"
"""

""" OG call
#Reading all the input files:
cns = "/data/Compass/dev/valdezkm/cnvkit/TSO_vs_ExomeList/CP04210_T1D_E/CP04210_T1D_E.bwa.final.call.cns"
#input_file = "CP04210_T1D_E.call.bed"
bed_file_input= "/data/Compass/Ref/hg19/design/NextraExome_v1.target.hg19.merged.bed"
gene_list = "/data/Compass/Ref/hg19/annotation/geneLists/combinedList_05312018"

#Reading all the temp files:
segment_temp1= "outputs/New_DefaultLowCov/CP04210_T1D_E_temp1.txt"
py_intersection1= "outputs/New_DefaultLowCov/CP04210_T1D_E_temp2.txt"
result_temp1= "outputs/New_DefaultLowCov/CP04210_T1D_E_temp3.txt"
gene_seperated_file= "outputs/New_DefaultLowCov/CP04210_T1D_E_temp4.txt"
py_intersection2 = "outputs/New_DefaultLowCov/CP04210_T1D_E_temp5.txt"
final_temp= "outputs/New_DefaultLowCov/CP04210_T1D_E_temp6.txt"
output_txt = "outputs/New_DefaultLowCov/CP04210_T1D_E.txt"
output_gl = "outputs/New_DefaultLowCov/CP04210_T1D_E_genelevel.txt"
#cns = "CP04210_T1D_E.call.cns"
cns_temp1 = "cns.tmp"
"""

choose_cn = lambda x: np.round(np.median(x[~np.isnan(x)]), 0)
gene_aliases = {
    "AFDN": "MLLT4",
    "NSD2": "WHSC1",
    "NSD3": "WHSC1L1",
    "SEPTIN9": "SEPT9",
    "TENT5C": "FAM46C",
    "H3-3A": "H3F3A",
    "H3C2": "HIST1H3B",
}
reverse_aliases = {v: k for k, v in gene_aliases.items()}
gene_aliases["TERC"] = np.nan


def gene_agg(data):
    aggs_res = (
        data.sort_values(by=["chromosome", "start", "end"])
        .groupby("gene", sort=False, as_index=False)
        .agg(
            {
                "chromosome": "first",
                "start": "first",
                "end": "last",
                "cn": "median",
                "cn1": choose_cn,
                "cn2": choose_cn,
            }
        )
    )
    # aggs_res['CNt']= aggs_res['CNt'].round(0).astype(int)
    # aggs_res['A'] = aggs_res['A'].round(0).astype(int)
    # aggs_res['B'] = aggs_res['B'].round(0).astype(int)
    # if A + B doesn't = CNt, then let's make A and B NAs
    for i in aggs_res.index:
        if (
            not aggs_res.loc[i, "cn"] == aggs_res.loc[i, "cn1"] + aggs_res.loc[i, "cn2"]
        ):  # or not aggs_res.loc[i,'A'] >= aggs_res.loc[i,'B']:
            aggs_res.loc[i, "cn1"], aggs_res.loc[i, "cn2"] = ("NA", "NA")
    return aggs_res


"""

# Put CNS in py_intersection1 format.

with open(input_file, 'r') as f, open(py_intersection1,'w') as fo:
    for line in f.readlines():
        if line[:6] == 'chromo':
            pass
        else:
            fields = line.split('\t')
            fields[3] = re.sub(r'(\d),([a-zA-Z])', r'\1;\2', fields[3])
            for gene in fields[3].split(';'):
                fo.write('\t'.join([fields[0],fields[1],fields[2],gene,'\n']))

"""


def __main__():
    with open(cns, "r") as f, open(cns_temp1, "w") as t:
        for line in f.readlines():
            fields = line.split("\t")
            if fields[0] == "chromosome":
                cns_header = [i.rstrip() for i in fields]
                cns_header_line = line
                continue
            t.write("\t".join(fields[0:3] + fields[4:]))

    # Removing quotes from the input file
    # with open(input_file, 'r') as f, open(segment_temp1, 'w') as fo:
    with open(cns_temp1, "r") as f, open(segment_temp1, "w") as fo:
        headerline = f.readline()
        headerline = headerline.replace('"', "").replace("'", "")
        # fo.write("#" + headerline + "\n")
        fo.write(headerline + "\n")
        fo.write(headerline + "\n")
        for line in f:
            fo.write(line.replace('"', "").replace("'", ""))

    # bedtools intersection using pybedtools of the bed_interval file and input file without quotations (segment_temp1)
    bed_file = pybedtools.BedTool(bed_file_input)
    segment_file = pybedtools.BedTool(segment_temp1)
    # CDKN2B_bed = bed_file.filter(lambda x: "CDKN2B" in x.name)

    result_temp = bed_file.intersect(segment_file, wa=True)
    result = result_temp.saveas(py_intersection1)

    # Replacing all the not found genes and extracting the 1st four columns
    with open(py_intersection1, "r") as f, open(result_temp1, "w") as final:
        final.write(
            "Chrom" + "\t" + "Start.pos" + "\t" + "End.pos" + "\t" + "Gene" + "\n"
        )
        list_random = []
        for line in f:
            line_without = line.replace("___", "\t")
            list_random = line_without.split()
            if not "NOTFOUND" in line_without:
                final.write("\t".join(list_random[0:4]) + "\n")

    # Spliting all the multiple genes value to different rows
    df = pd.read_csv(result_temp1, delimiter="\s+", index_col=False)
    df = df.drop("Gene", axis=1).join(
        df["Gene"]
        .str.split(",", expand=True)
        .stack()
        .reset_index(level=1, drop=True)
        .rename("Gene")
    )
    df.to_csv(gene_seperated_file, header=None, index=None, sep="\t")

    # Bedtools interesection with the newgene separated file and the segment file
    # temp1_file = pybedtools.BedTool(gene_seperated_file)

    # result_intersection = segment_file.intersect(temp1_file,wb = True)
    # result = result_intersection.saveas(py_intersection2)

    # Mimic with cns:
    # first remove genes from cns
    # cns_df = pd.read_csv(cns, sep = '\t')
    # cns_df[[i for i in cns_df.columns if not i == 'gene']].to_csv(,sep = '\t')
    # Bedtools interesection with the newgene separated file and the cns file
    temp1_file = pybedtools.BedTool(gene_seperated_file)
    cns_bed = pybedtools.BedTool(cns_temp1)

    cns_intersect = cns_bed.intersect(temp1_file, wb=True)
    cns_results = cns_intersect.saveas(py_intersection2)

    # Final temp file with all the bed and list details
    line_final = []
    genes = []
    with open(py_intersection2, "r") as f, open(final_temp, "w") as fo:
        # fo.write('#chromosome\tstart.pos\tend.pos\tGene\tBf\tN.BAF\tsd.BAF\tdepth.ratio\tN.ratio\tsd.ratio\tCNt\tA\tB\tLPP\tGene\tSource;Source\t#Lists\n')
        for line in f:
            line = line.strip()
            list_final = line.split("\t")
            # Maintain NAs when BAF, A and B are blank
            for i in range(len(list_final)):
                if list_final[i] == "":
                    list_final[i] = "NA"
            gene_index = list_final[-1:]
            genes.append(gene_index)
            first_index = list_final[0:3]
            middle_index = list_final[3:-4]
            #### Fill in the NAs in As and Bs
            newline = (
                "\t".join(first_index)
                + "\t"
                + "\t".join(gene_index)
                + "\t"
                + "\t".join(middle_index)
                + "\n"
            )
            fo.write(newline)

    # Check if the temporary file is empty
    if os.stat(final_temp).st_size != 0:
        # Reading the combined gene list text file and converting into a dictionary
        # gene_list = "/data/Compass/Ref/hg19/annotation/geneLists/combinedList_05312018"
        ann_dict = {}
        ann_array = []
        with open(gene_list) as ann_read:
            for line in ann_read:
                ann_array = line.split()
                key_gene = ann_array[0]
                ann_dict[key_gene] = ann_array[1], ann_array[2]
        # Annotation of the genes
        temp_array = []
        OutputFile = open(output_txt, "w")
        # OutputFile.write('#chromosome\tstart.pos\tend.pos\tGene\tBf\tN.BAF\tsd.BAF\tdepth.ratio\tN.ratio\tsd.ratio\tCNt\tA\tB\tLPP\tgene_annotation\tSource;Source\t#Lists\n')
        # OutputFile.write('#chromosome\tstart.pos\tend.pos\tGene\tlog2\tbaf\tci_hi\tci_lo\tCNt\tA\tB\tdepth\tprobes\tweight\tgene_annotation\tSource;Source\t#Lists\n')
        OutputFile.write(
            f"{cns_header_line.rstrip()}\tgene_annotation\tSource;Source\t#Lists\n"
        )
        with open(final_temp, "r") as temp1_read:
            # next(temp1_read)
            for line in temp1_read:
                line = line.strip()
                temp_array = line.split("\t")
                gene_value = temp_array[3]
                if gene_value in ann_dict.keys():
                    list_values = list(ann_dict[gene_value])
                    if gene_value in reverse_aliases.keys():
                        temp_array[3] = reverse_aliases[gene_value]
                    OutputFile.write(
                        "\t".join(temp_array)
                        + "\t"
                        + gene_value
                        + "\t"
                        + "\t".join(list_values)
                        + "\n"
                    )
                else:
                    OutputFile.write(
                        "\t".join(temp_array)
                        + "\t"
                        + "-"
                        + "\t"
                        + "-"
                        + "\t"
                        + "-"
                        + "\n"
                    )
        OutputFile.close()
        # aggregating at the gene level and taking the median for CNV
        df_read = pd.read_csv(output_txt, delimiter="\t")
        df_1 = gene_agg(df_read)
        df_final = df_1[["chromosome", "start", "end", "gene", "cn", "cn1", "cn2"]]
        df_final["cn"] = np.floor(
            pd.to_numeric(df_final["cn"], errors="coerce")
        ).astype("Int64")
        df_final["cn1"] = np.floor(
            pd.to_numeric(df_final["cn1"], errors="coerce")
        ).astype("Int64")
        df_final["cn2"] = np.floor(
            pd.to_numeric(df_final["cn2"], errors="coerce")
        ).astype("Int64")
        df_final.loc
        df_final.to_csv(output_gl, sep="\t", header=True, index=False)
    else:
        print("The temp file is empty... Probably due to an empty .cns")
        with open(output_gl, "w") as o:
            pass


outfile1 = cns.split("/")[-1]
with contextlib.redirect_stdout(io.StringIO()):
    __main__()


print(f"{sys.argv[0]} complete for {outfile1}")
