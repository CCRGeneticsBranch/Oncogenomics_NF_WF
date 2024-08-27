There are two ways to generate a samplesheet for the pipeline depending on the use case. If the goal is to process patient fastq files through the pipeline to use the results for secondary analysis, you can follow these steps to [build your own samplesheet](#build-your-own-samplesheet).
Along with processing the data, if you want to visualize the results on [ClinOmics data portal](https://oncogenomics.ccr.cancer.gov/production/public/) then follow steps to [build samplesheet from mastersheet](#build-samplesheet-from-mastersheet). This is highly recommended.

## Build samplesheet from mastersheet

For khanlab purposes, pipeline is always launched using the information in the mastersheets on biowulf under /data/khanlab space. The script [samplesheet_builder.py](https://github.com/CCRGeneticsBranch/Oncogenomics_NF_WF/blob/feature/casename/samplesheet_builder.py) queries the mastersheets to build a samplesheet for the pipeline. A copy of this script is available in the pipeline git repo. This script takes two inputs `PatientID` and `casename`. By default, it queries all mastersheets found in the `/data/khanlab/projects/DATA/Sequencing_Tracking_Master` directory and uses `/data/khanlab/projects/DATA` as the default input directory.

```
Usage: python ./samplesheet_builder.py <patient_id> <case_name>
Default Samplesheet Directory: /data/khanlab/projects/DATA/Sequencing_Tracking_Master
Default Input Directory: /data/khanlab/projects/DATA
To use custom directories, modify the script:
   - Change 'DEFAULT_SAMPLESHEET_DIR' to your samplesheet directory path
   - Change 'DEFAULT_INPUT_DIR' to your input directory path
```

` python ./samplesheet_builder.py Test_Patient casename` will output a file Test_Patient_casename.csv in the same folder.

## Build your own samplesheet

Alternately, we can build the pipeline input csv file using these columns.

| Column name     | Notes                                                                                             | Example                                                                            |
| --------------- | ------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| sample          | Patient name                                                                                      | NCI-Test1                                                                          |
| library         | Name of the sample library                                                                        | Test1_T1D_E                                                                        |
| read1           | Full path to the read1                                                                            | /data/khanlab/DATA/Sample_Test1_T1D_E/Sample_Test1_T1D_E.R1.fastq.gz               |
| read2           | Full path to the read2                                                                            | /data/khanlab/DATA/Sample_Test1_T1D_E/Sample_Test1_T1D_E.R2.fastq.gz               |
| sample_captures | Name of the capture kit used                                                                      | List of supported capture kits are [here](index.md#sequencing-capture-kits)        |
| Matched_RNA     | Matched RNA library for the tumor library. This includes cell_line_RNA and tumor_RNA              | Test1_T1R_T                                                                        |
| Matched_normal  | Matched exome normal library for the tumor library. This includes panel, blood DNA, cell_line_DNA | Test1_N1D_E                                                                        |
| Diagnosis       | Diagnosis of the patient                                                                          | Glioma                                                                             |
| casename        | Casename for the patient                                                                          | NCI-Test1                                                                          |
| type            | Data type                                                                                         | example: tumor_RNA, tumor_DNA, normal_DNA, blood_DNA, cell_line_DNA, cell_line_RNA |
| FCID            | Flowcell ID                                                                                       | ACJ678349                                                                          |

### example samplesheet:

sample,library,read1,read2,sample_captures,Diagnosis,Matched_RNA,Matched_normal,casename,type,FCID,Project
Test8,Test5_T1D_E,/data/khanlab/projects/fastq/Test5_T1D_E_R1.fastq.gz,/data/khanlab/projects/fastq/Test5_T1D_E_R2.fastq.gz,clin.ex.v1,Osteosarcoma,,Test8_N2D_E,NFtest0523,tumor_DNA,AWXYNH2,Test
Test8,Test8_N2D_E,/data/khanlab/projects/fastq/Test8_N2D_E_R1.fastq.gz,/data/khanlab/projects/fastq/Test8_N2D_E_R2.fastq.gz,clin.ex.v1,Osteosarcoma,,,NFtest0523,normal_DNA,AWXYNH2,Test
