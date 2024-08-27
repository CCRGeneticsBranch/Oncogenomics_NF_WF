# Running the workflow.

To initiate the workflow, execute the `launch.sh` script in a interactive node. This script acts as a wrapper, spawning and submitting jobs to the SLURM scheduling system. It extracts the `PatientID` and `Casename` information from the [samplesheet](samplesheet.md), creates the corresponding output directory for each case, and organizes all results and log files under `[output_directory]/PatientID/Casename/`.

```
Usage: ./launch.sh <samplesheet_with_full_path> [output_directory] [genome]
This script requires at least one positional argument:
1. Full path to samplesheet csv file
Optional arguments:
2. Path to results directory (default: /data/khanlab/projects/processed_DATA)
3. Genome name. Accepted values are hg19 and mm39 (default: hg19)

```

# Workflow log

When the workflow is launched, it will produce a log file that provides information about the pipeline execution, including the command line used, the version of Nextflow, the input folder path, the results directory, and the work directory.

The log file will be generated under /Resultsdir/PatientID/casename directory. This is the logfile naming convention. `patientid_casename_jobid_datelaunched-timelaunched.out`

`Example: NCI0439_TestTNR_34233573_20240822-183525.out`

### Sample log file:

```
[+] Loading singularity  4.0.3  on cn4280
[+] Loading java 17.0.3.1  ...
[+] Loading nextflow  23.10.0
[+] Loading Graphviz v 2.46.1  ...
NXF_HOME=/data/khanlab/projects/Nextflow_dev/dev/NCI0439/TestTNR/.nextflow
nextflow run -c /vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/nextflow.config -profile biowulf_test_run_slurm --logdir /data/khanlab/projects/Nextflow_dev/dev/NCI0439/TestTNR/log /vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/main.nf -resume -with-trace -with-timeline -with-dag
ESC[33mNextflow 24.04.4 is available - Please consider updating your version to itESC(BESC[m
N E X T F L O W  ~  version 23.10.0
Launching `/vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/main.nf` [shrivelled_stallman] DSL2 - revision: b8c8cd72d2
E X O M E - R N A S E Q - N F   P I P E L I N E
===================================
NF version   : 23.10.0
runName      : shrivelled_stallman
username     : gangalapudiv2
configs      : [/vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/nextflow.config, /vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/nextflow.config]
profile      : biowulf_test_run_slurm
cmd line     : nextflow run -c /vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/nextflow.config -profile biowulf_test_run_slurm --logdir /data/khanlab/projects/Nextflow_dev/dev/NCI0439/TestTNR/log /vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/main.nf -resume -with-trace -with-timeline -with-dag
start time   : 2024-08-20T12:54:52.989621486-04:00
projectDir   : /vf/users/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF
launchDir    : /vf/users/khanlab/projects/Nextflow_dev/dev/NCI0439/TestTNR
workdDir     : /vf/users/khanlab/projects/Nextflow_dev/dev/NCI0439/TestTNR/work
homeDir      : /home/gangalapudiv2

[-        ] process > Tumor_Normal_RNAseq_WF:Comm... -
[-        ] process > Tumor_Normal_RNAseq_WF:Comm... -
[-        ] process > Tumor_Normal_RNAseq_WF:Comm... -
[-        ] process > Tumor_Normal_RNAseq_WF:Comm... -
[-        ] process > Tumor_Normal_RNAseq_WF:Comm... -

```

# Workflow Resources

$Script_home ---> Path to the code repository (/data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF).

$Script_home/nextflow.config ---> All the pipeline config resources are called from this file. We are using `biowulf_test_run_slurm` profile to run the samples on biowulf.

$Script_home/config ---> Cluster config, singularity config, docker config and params config are located here.

$Script_home/main.nf ---> This script calls the workflows and subworkflows in the pipeline.

# Run Debugging

Run-related Files: All files associated with the run will be generated in the `[output_directory]/PatientID/Casename` directory.

Log Files: The log file, ending in `*.out`, will be located in `[output_directory]/PatientID/Casename`. This file contains details on any errors encountered during the run and provides the path to the work directory for troubleshooting.

Work directory: Inside `[output_directory]/PatientID/Casename/work`, you will find subfolders containing `.command.sh`, `.command.log`, and `.command.err` files for each process. These files offer detailed information about the execution of each command and any issues that may have occurred.
