## AWS_POC_MVP_NF

We started the [AWS_MVP_HPC](https://github.com/CCRGeneticsBranch/AWS_MVP_HPC) project to kick start transfer of CCR Genetics Branch (GB) NGS data analyses pipelines to the cloud (specifically AWS AGC.) Since, AWS(Cloud One) [AGC](https://aws.github.io/amazon-genomics-cli/docs/reference/agc/) supports multiple pipelining frameworks we decided to test out two:
 - Nextflow
 - Snakemake

This code repo represents the Nextflow implementation of the POC (and MVP) of the GB's transcriptomic pipeline currently running on [BIOWULF](https://hpc.nih.gov/) using Snakemake.

### Flowchart

https://i.imgur.com/NvfwmRD.png

This chart represents the entire MVP plan with the POC being highlighted in yellow.

> <ins>Disclaimer</ins>: There are two parallel efforts underway for the _"cloudification"_ of GB, namely, 
> 
> a. Moving the database management and its web-interface to AWS also referred to as **"the Database component"** and 
> 
> b. Orchestrating NGS analysis workflows on AWS which is also referred to as **"the HPC compotent"**. This repository solely focuses on _the HPC component_.

For complete documentation please go [here](https://CCRGeneticsBranch.github.io/AWS_MVP_HPC/).

Please send your comments/questions/suggestions to [Vishal Koparde](https://github.com/kopardev) via [email](mailto:vishal.koparde@nih.gov).

