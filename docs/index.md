## Overview of the workflow
This is Nextflow based RNAseq workflow that performs read trimming, paired end STAR alignment, quantification of transcript abundances, Fusion calling, GATK based 
variant analysis, Annotation using Annovar and comprehensive QC. 
![RNAseq_workflow](DAG_rnaseq.png)

# Prerequisites

This workflow requires modules Nextflow, Singularity and graphviz


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
directory. This workflow expects all the fastq files in same directory.
3. Edit the resources config file `AWS_POC_Nextflow/config/Run_upto_quants_cluster.config` and update the cluster requirements for the processes.
4. Singularity config is set up `AWS_POC_Nextflow/config/biowulf_singularity.config` please update the `runOptions` as per your cluster environment.
5. `nextflow.config` has a profile `Run_upto_quants` that use the config files listed in steps 2,3,4. 
6. Finally update the path for the variable `export OUTDIR` in the script `Run_upto_quants.sh`. 


# Running the workflow.

 The workflow can be launched using the Run_upto_quants.sh script. 
 ```
This script takes two inputs
Provide Tag argument - This will add a tag to your resultsdir.
Provide Run_upto_counts_only argument - It takes value true/false
example: sh nf.sh projectname true #this will run upto RSEM and generates counts 
example: sh nf.sh projectname false #this will run the complete pipeline 

 ```
The following command should kick start the pipeline and display the status.

`sh Run_upto_quants.sh Neuroblastoma_MAQC-III_SEQC true`

On launching the workflow,it prints out a similar log with information on the command line used to run the pipeline, Nextflow version, input folder path, results 
directory and work directory

```
+] Loading singularity  3.10.5  on cn4271 
[+] Loading java 17.0.2  ... 
[+] Loading nextflow 22.10.2
[+] Loading Graphviz v 2.40.1  ... 
NXF_HOME=/data/khanlab3/kopardevn/AWS_MVP_test/results.vg1.2/.nextflow
nextflow run -c /vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/nextflow.config -profile Run_upto_quants 
/vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/main.nf -resume --run_upto_counts true -with-trace -with-timeline 
/data/khanlab3/kopardevn/AWS_MVP_test/results.vg1.2/timeline.html
N E X T F L O W  ~  version 22.10.4
Launching `/vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/main.nf` [crazy_tuckerman] DSL2 - revision: 3e49210bb9
R N A S E Q - N F   P I P E L I N E  
===================================
NF version   : 22.10.4
runName      : crazy_tuckerman
username     : gangalapudiv2
configs      : [/vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/nextflow.config, /vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/nextflow.config]
profile      : Run_upto_quants
cmd line     : nextflow run -c /vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/nextflow.config -profile Run_upto_quants 
/vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/main.nf -resume --run_upto_counts true -with-trace -with-timeline 
/data/khanlab3/kopardevn/AWS_MVP_test/results.vg1.2/timeline.html
start time   : 2023-02-09T13:39:52.890203954-05:00
projectDir   : /vf/users/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow
launchDir    : /vf/users/khanlab3/kopardevn/AWS_MVP_test/results.vg1.2
workdDir     : /data/khanlab3/kopardevn/AWS_MVP_test/work.vg1.2
homeDir      : /home/gangalapudiv2
reads        : /data/khanlab/projects/fastq/*_{R1,R2}.fastq.gz

[b5/33c3ef] process > Cutadapt (Test3_R_T)       [100%] 3 of 3, cached: 3 ✔
[70/1fd193] process > Fastqc (Test3_R_T)         [100%] 3 of 3, cached: 3 ✔

```


