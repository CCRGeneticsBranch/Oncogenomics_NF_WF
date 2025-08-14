## Usage


This page describes how to run the workflow on **Biowulf**.  
For installation and prerequisites, see the [Installation](index.md) page.

---

## Running on Biowulf

On Biowulf, the pipeline is hosted in the Khanlab space.  
A helper script, `launch.sh`, is provided to simplify execution. This script wraps the `nextflow run` command with Biowulf-specific settings and automatically configures output directories, logging, and profiles.

### Step 1 — Copy the launch script

??? "launch.sh"

    ```
    #!/bin/bash

    # Set default values
    DEFAULT_OUTDIR="/data/khanlab/projects/processed_DATA"
    DEFAULT_GENOME="hg19"
    DEFAULT_PLATFORM="biowulf"

    # Check if the required samplesheet argument is provided
    if [[ "$#" -lt 1 ]]; then
        echo "Usage: $0 <samplesheet_with_full_path> [output_directory] [genome]"
        echo "This script requires at least one positional argument:"
        echo "1. Path to samplesheet"
        echo "Optional arguments:"
        echo "2. Path to results directory (default: $DEFAULT_OUTDIR)"
        echo "3. Genome name. Accepted values are hg19 and mm39 (default: $DEFAULT_GENOME)"
        exit 1
    fi

    # Assign the first argument (samplesheet) to a variable
    export SAMPLESHEET=$1

    # Assign the second argument (output_directory) if provided, otherwise use the default
    export OUTDIR=${2:-$DEFAULT_OUTDIR}

    # Assign the third argument (genome) if provided, otherwise use the default
    export GENOME=${3:-$DEFAULT_GENOME}

    export PLATFORM=${4:-$DEFAULT_PLATFORM}

    # Ensure the genome name is valid
    if [[ "$GENOME" != "hg19" && "$GENOME" != "mm39" ]]; then
        echo "Invalid genome specified. Accepted values are hg19 and mm39."
        exit 1
    fi

    WF_HOME="/data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF"
    CONFIG_FILE="$WF_HOME/biowulf_nextflow.config"

    #export PATIENT=$(awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="sample") s=i} NR>1 {print $s}' "$SAMPLESHEET" | sort | uniq)
    #export CASENAME=$(awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="casename") c=i} NR>1 {print $c}' "$SAMPLESHEET" | sort | uniq)
    export PATIENT=$(python3 -c 'import csv,sys; r=csv.DictReader(open(sys.argv[1])); print("\n".join(sorted(set(row["sample"] for row in r))))' "$SAMPLESHEET")
    export CASENAME=$(python3 -c 'import csv,sys; r=csv.DictReader(open(sys.argv[1])); print("\n".join(sorted(set(row["casename"] for row in r))))' "$SAMPLESHEET")

    export RESULTSDIR="$OUTDIR/$PATIENT/$CASENAME"

    mkdir -p "$RESULTSDIR"

    export LOG="$RESULTSDIR/log"

    mkdir -p "$LOG"

    export NXF_HOME="$RESULTSDIR/.nextflow"

    if [[ -z "$PATIENT" || -z "$CASENAME" ]]; then
        echo "Error: Could not extract PATIENT or CASENAME from the samplesheet."
        exit 1
    fi


    cd $RESULTSDIR

    if [[ "$GENOME" == "hg19" ]]; then
        PROFILE="biowulf_test_run_slurm"
    elif [[ "$GENOME" == "mm39" ]]; then
        PROFILE="biowulf_mouse_RNA_slurm"
    else
        echo "Unknown genome: $GENOME"
        exit 1
    fi

    logname=$(basename "$SAMPLESHEET" .csv)
    timestamp=$(date +"%Y%m%d-%H%M%S")
    sbatch <<EOT
    #!/bin/bash
    #SBATCH --job-name="$logname"
    #SBATCH --output="$OUTDIR/$PATIENT/$CASENAME/${logname}_%A_${timestamp}.out"
    #SBATCH --cpus-per-task=2
    #SBATCH --mem=05g
    #SBATCH --time=08-00:00:00

    module load nextflow/23.10.0 singularity graphviz
    nextflow run -c $CONFIG_FILE -profile $PROFILE --logdir $LOG $WF_HOME/main.nf -resume --samplesheet $SAMPLESHEET --resultsdir $OUTDIR --genome_v $GENOME --platform $PLATFORM
    exit 0
    EOT

    ```

Download or copy the script into your Biowulf working directory and make it executable:

```
cp /data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/launch.sh .
chmod +x launch.sh
```

### Step 2 — View usage instructions
Run the script without arguments to display the usage help:

```
./launch.sh

Usage: ./launch.sh <samplesheet_with_full_path> [output_directory] [genome]
This script requires at least one positional argument:
1. Path to samplesheet
Optional arguments:
2. Path to results directory (default: /data/khanlab/projects/processed_DATA)
3. Genome name. Accepted values are hg19 and mm39 (default: hg19)
```

Usage Examples:

1. Run with defaults (output to default path, genome hg19):

    ```./launch.sh /path/to/samplesheet.csv```

2. Run with custom output directory:

    ```./launch.sh /path/to/samplesheet.csv /custom/output/path```

3. Run with custom output directory and genome:

    ```./launch.sh /path/to/samplesheet.csv /custom/output/path mm39```

!!! tip "mm39"
    If mm39 is selected for the genome, the pipeline will run *Mapping & Gene expression* only for mouse data.


## Preparing the Samplesheet

You can prepare a samplesheet in two ways, depending on your use case:

1. **Build from a Master Sheet** – Recommended if you also plan to visualize results in the [ClinOmics Data Portal](https://oncogenomics.ccr.cancer.gov/production/public/).  
2. **Build Manually** – For standalone processing without portal integration.

### 1. Build from a Master Sheet (Recommended)

For Khanlab workflows, master sheets are stored on Biowulf in:

```/data/khanlab/projects/DATA/Sequencing_Tracking_Master```

Use the script:

[samplesheet_builder.py](https://github.com/CCRGeneticsBranch/Oncogenomics_NF_WF/blob/feature/casename/samplesheet_builder.py)
```
python ./samplesheet_builder.py <PatientID> <CaseName>
```

Default Master Sheet Directory: /data/khanlab/projects/DATA/Sequencing_Tracking_Master

Default Input Directory: /data/khanlab/projects/DATA

To use custom directories, edit DEFAULT_SAMPLESHEET_DIR and DEFAULT_INPUT_DIR in the script.

**Error Handling:**

1. **Invalid FASTQ Paths** – If `read1` or `read2` paths are invalid, the script will print an error prompting you to verify input paths.
2. **Unknown Capture Kit** – If the capture kit entered is not listed in [Sequencing Capture Kits](index.md#sequencing-capture-kits), the script will assign `unknown` as the kit name.

### 2. Build Your Own Samplesheet

Create a CSV with these column information:

| Column          | Description                                        | Example                                                                 |
|-----------------|----------------------------------------------------|-------------------------------------------------------------------------|
| sample          | Patient name                                       | NCI-Test1                                                               |
| library         | Sample library name                                | Test1_T1D_E                                                             |
| read1           | Path to R1 FASTQ                                   | /data/khanlab/DATA/Sample_Test1_T1D_E/Sample_Test1_T1D_E.R1.fastq.gz    |
| read2           | Path to R2 FASTQ                                   | /data/khanlab/DATA/Sample_Test1_T1D_E/Sample_Test1_T1D_E.R2.fastq.gz    |
| sample_captures | Capture kit name                                   | [Sequencing Capture Kits](index.md#sequencing-capture-kits)                                                               |
| Matched_RNA     | Matched RNA library (optional)                     | Test1_T1R_T                                                             |
| Matched_normal  | Matched exome normal library (optional)            | Test1_N1D_E                                                             |
| Diagnosis       | Patient diagnosis                                  | Glioma                                                                  |
| casename        | Case name                                          | NCI-Test1                                                               |
| type            | Data type                                          | tumor_RNA, tumor_DNA, normal_DNA, etc.                                  |
| FCID            | Flowcell ID (optional)                             | ACJ678349                                                               |
| Genome          | Genome                                             | hg19 or mm39                                                            |


**Example Input csv:**

| sample | library      | read1                                                                                       | read2                                                                                       | sample_captures | Diagnosis    | Matched_RNA | Matched_normal | casename   | type       | FCID     | Genome |
|--------|-------------|---------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|-----------------|--------------|-------------|----------------|------------|------------|----------|---------|
| Test8  | Test5_T1D_E | /data/khanlab/projects/fastq/Test5_T1D_E_R1.fastq.gz                                         | /data/khanlab/projects/fastq/Test5_T1D_E_R2.fastq.gz                                         | clin.ex.v1      | Osteosarcoma |             | Test8_N2D_E    | NFtest0523 | tumor_DNA | AWXYNH2  | hg19    |
| Test8  | Test8_N2D_E | /data/khanlab/projects/fastq/Test8_N2D_E_R1.fastq.gz                                         | /data/khanlab/projects/fastq/Test8_N2D_E_R2.fastq.gz                                         | clin.ex.v1      | Osteosarcoma |             |                | NFtest0523 | normal_DNA| AWXYNH2  | hg19    |

