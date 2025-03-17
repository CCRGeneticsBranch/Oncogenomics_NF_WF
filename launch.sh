#!/bin/bash

# Set default values
DEFAULT_OUTDIR="/data/khanlab/projects/processed_DATA"
DEFAULT_GENOME="hg19"

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

# Ensure the genome name is valid
if [[ "$GENOME" != "hg19" && "$GENOME" != "mm39" ]]; then
    echo "Invalid genome specified. Accepted values are hg19 and mm39."
    exit 1
fi

WF_HOME="/data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF"
CONFIG_FILE="$WF_HOME/nextflow.config"

export PATIENT=$(awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="sample") s=i} NR>1 {print $s}' "$SAMPLESHEET" | sort | uniq)
export CASENAME=$(awk -F',' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="casename") c=i} NR>1 {print $c}' "$SAMPLESHEET" | sort | uniq)
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
nextflow run $WF_HOME/main.nf -resume --samplesheet $SAMPLESHEET --resultsdir $OUTDIR --genome_v $GENOME -profile $PROFILE --logdir $LOG
exit 0
EOT
