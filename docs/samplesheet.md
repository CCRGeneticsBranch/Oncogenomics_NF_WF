There are two ways to generate a samplesheet for the pipeline depending on the use case. If the goal is to process patient fastq files through the pipeline to use the results for secondary analysis, you can follow these steps to [build your own samplesheet](#build-your-own-samplesheet).
Along with processing the data, if you want to visualize the results on [ClinOmics data portal](https://oncogenomics.ccr.cancer.gov/production/public/) then follow steps to [build samplesheet from mastersheet](#build-samplesheet-from-mastersheet). This is highly recommended.

## Build your own samplesheet

The pipeline requires a CSV file as input for each patient. The samplesheet must include the following mandatory columns:

| Column name     | Notes                                              | Example                                                                     |
| --------------- | -------------------------------------------------- | --------------------------------------------------------------------------- |
| sample          | Patient name                                       | NCI-Test1                                                                   |
| library         | Name of the sample library                         | Test1_T1D_E                                                                 |
| read1           | Full path to the read1                             | /data/khanlab/DATA/Sample_Test1_T1D_E/Sample_Test1_T1D_E.R1.fastq.gz        |
| read2           | Full path to the read2                             | /data/khanlab/DATA/Sample_Test1_T1D_E/Sample_Test1_T1D_E.R1.fastq.gz        |
| sample_captures | Name of the capture kit used                       | List of supported capture kits are [here](index.md#sequencing-capture-kits) |
| Matched_RNA     | Matched RNA library for the tumor library          | Test1_T1R_T                                                                 |
| Matched_normal  | Matched exome normal library for the tumor library | Test1_N1D_E                                                                 |
| Diagnosis       | Diagnosis of the patient                           | Glioma                                                                      |
| casename        | Casename for the patient                           | NCI-Test1                                                                   |
| type            | Data type                                          | example: tumor_RNA, tumor_DNA, normal_DNA, blood_DNA                        |
| FCID            | Flowcell ID                                        | ACJ678349                                                                   |

## Build samplesheet from mastersheet
