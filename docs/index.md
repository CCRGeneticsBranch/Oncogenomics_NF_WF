## Overview of the workflow
This is Nextflow based RNAseq workflow that performs read trimming, paired end STAR alignement, quantification of transcript abundances, Fusion calling, GATK based 
variant analysis, Annotation using Annovar and comprehensive QC. 
![RNAseq_workflow](/Applications/Projects/DAG_rnaseq.png)

# Installation
Please clone this repository to your local filesystem using the following command:

```
        git clone https://github.com/CCRGeneticsBranch/AWS_POC_Nextflow.git
        cd AWS_POC_Nextflow/ 
        git checkout feature/mvp
```

# Setting up the workflow to run upto quants. 
1. Download the References folder from [box](https://nih.app.box.com/folder/193831680410) to a accessable location and untar `index-STAR_2.7.9a.tar.gz` and 
`rsem_1.3.2.tar.gz` file.
2. Edit the references config file  `AWS_POC_Nextflow/config/Run_upto_quants_references.config` and update the file paths of the references and path to the input 
directory. This worflow expects all the fastq files in same directory.
3. Edit the resources config file `AWS_POC_Nextflow/config/Run_upto_quants_cluster.config` and update the cluster requirements for the processes.
4. Singulairty config is set up `AWS_POC_Nextflow/config/biowulf_singularity.config` please update the `runOptions` as per your cluster environment.
5. `nextflow.config` has a profile `Run_upto_quants` to uses the config files listed in steps 2,3,4. 


# Running the workflow.

 The workflow can be launched using the Run_upto_quants.sh script. 
 ```
This script takes two inputs
Provide Tag argument - This will add a tag to your resultsdir.
Provide Run_upto_counts_only argument - It takes value true/false
example: sh nf.sh projectname true #this will run upto RSEM and generates counts 
example: sh nf.sh projectname fasle #this will run the complete pipeline 

 ```
The following command should kick start the pipeline and display the status.

`sh Run_upto_quants.sh Neuroblastoma_MAQC-III_SEQC true`

