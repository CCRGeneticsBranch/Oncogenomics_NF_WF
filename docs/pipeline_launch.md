# Running the workflow.

This page explains what gets written where, how to monitor runs on Biowulf, how to resume or recover failed tasks, and what to collect when reporting errors.



### Results & Directory Layout:

${OUTDIR} is the results root — either the default location or the custom path you provide when [launching](usage.md) the workflow.

| Files | Path |
|---|---|
| Results root | `${OUTDIR}/${PATIENT}/${CASENAME}/` |
| NXF state | `${OUTDIR}/${PATIENT}/${CASENAME}/.nextflow/` |
| Run logs (stdout) | `${OUTDIR}/${PATIENT}/${CASENAME}/*.out` |
| Nextflow reports (timeline, trace, manifest) | `${OUTDIR}/${PATIENT}/${CASENAME}/log/` |
| Nextflow workdir | `${OUTDIR}/${PATIENT}/${CASENAME}/work/` |


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


### Debugging: Where to look

#### 1) Start with the run’s stdout
Check the Biowulf/Slurm stdout files for a quick overview and the **exact work directory** of the failed process:

- `${OUTDIR}/${PATIENT}/${CASENAME}/*.out`

> These files usually contain the error summary and a line like `work/<hash>` pointing to the failing task’s work dir.

---

#### 2) Inspect per-process work-directory logs
For any failed task, open the work directory reported above and examine:

```
    work/<hash>/.command.sh
    work/<hash>/.command.out
    work/<hash>/.command.err
    work/<hash>/.command.log
```
> These files contain in-depth details of the process error (command, stdout, stderr, and Nextflow’s execution log).

---

### Reporting pipeline issues
If the failure appears to be a **code error** or a **pipeline glitch**, please open a new issue with logs [here](https://github.com/CCRGeneticsBranch/Oncogenomics_NF_WF/issues/new)