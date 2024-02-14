#!/usr/bin/env python
# script incorporated from COMPASS team
import pandas as pd
import sys
import re

opti = sys.argv[1]
HD = sys.argv[2]
out = sys.argv[3]


def checkData(dat):
    dataCheck = 0
    with open(dat) as f:
        file_content = f.read()
        if "A1" not in file_content and "HLA" not in file_content:
            dataCheck = 1
    return dataCheck


# optiType
optiCheck = checkData(opti)
if optiCheck == 0:
    opti = pd.read_csv(opti, sep="\t")
    opti = opti.drop(["Unnamed: 0", "Reads", "Objective"], axis=1)
    opti = opti.transpose()
    opti = opti.dropna()
    opti.columns = ["#Allele"]
    opti["#Allele"] = "HLA-" + opti["#Allele"]
    opti["OptiType"] = "Called"
    #    opti['OptiType'] = opti.index
    opti = opti.reset_index(drop=True)

# HLA-HD - has variable num of cols
HDcheck = checkData(HD)
if HDcheck == 0:
    HD = open(HD, "r")
    HD = [line.split("\t") for line in HD]
    HD = {line[0]: [item.strip() for item in line[1:]] for line in HD}
    HD = pd.DataFrame.from_dict(HD, orient="index")
    HD.loc[HD.iloc[:, 1] == "-", 1] = HD.iloc[:, 0]
    HD = HD.unstack().to_frame()
    HD.columns = ["#Allele"]
    HD = HD.dropna()
    HD = HD[HD["#Allele"].str.contains("HLA-A|HLA-B|HLA-C", regex=True)]
    HD["#Allele"] = HD["#Allele"].str.split(":").str[:2].str.join(":")
    HD["HLA-HD"] = "Called"
    HD = HD.reset_index(drop=True)

# if one fails, write other to output, otherwise merge and output
if optiCheck == 1 and HDcheck == 1:
    outF = open(out, "w")
    outF.write("#Allele\tOptiType\tHLA-HD")
elif optiCheck == 1 and HDcheck == 0:
    HD["OptiType"] = ""
    HD = HD[["#Allele", "OptiType", "HLA-HD"]]
    HD.to_csv(out, sep="\t", index=False)
elif optiCheck == 0 and HDcheck == 1:
    opti["HLA-HD"] = ""
    opti.to_csv(out, sep="\t", index=False)
else:
    # account for homozygous alleles so merge performs as expected
    opti["#Allele"][opti["#Allele"].duplicated()] = opti["#Allele"] + "_dup"
    HD["#Allele"][HD["#Allele"].duplicated()] = HD["#Allele"] + "_dup"
    consen = pd.merge(opti, HD, on="#Allele", how="outer")
    consen["#Allele"] = consen["#Allele"].str.replace("_dup", "", regex=False)
    consen.to_csv(out, sep="\t", index=False)
