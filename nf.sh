#!/bin/bash

if [[ "$#" -ne "5" ]]; then
    echo "This script requires five positional arguments:"
    echo "1. Path to samplesheet"
    echo "2. Path to results directory"
    echo "3. Patient name"
    echo "4. Case name"
    echo "5. Genome"
    exit 1
fi


set -e

SCRIPT_NAME="$0"
SCRIPT_DIRNAME=$(readlink -f $(dirname $0))
SCRIPT_BASENAME=$(basename $0)
WF_HOME=$SCRIPT_DIRNAME

CONFIG_FILE="$WF_HOME/nextflow.config"


# load singularity and nextflow modules
module load singularity nextflow/23.10.0 graphviz


export SAMPLESHEET=$1
export OUTDIR=$2
export PATIENT=$3
export CASENAME=$4
export GENOME=$5

# Check the value of the GENOME variable
if [[ "$GENOME" == "hg19" ]]; then
    PROFILE="biowulf_test_run_slurm"
elif [[ "$GENOME" == "mm39" ]]; then
    PROFILE="biowulf_mouse_RNA_slurm"
else
    echo "Unknown genome: $GENOME"
    exit 1
fi

export RESULTSDIR="$OUTDIR/$PATIENT/$CASENAME"

mkdir -p "$RESULTSDIR"

export WORKDIR="$RESULTSDIR/work"

export LOG="$RESULTSDIR/log"

mkdir -p "$LOG"

# set .nextflow dir ... dont want this to go to $HOME/.nextflow
export NXF_HOME="$RESULTSDIR/.nextflow"
# export SINGULARITY_BIND="/lscratch/$SLURM_JOB_ID"

printenv|grep NXF

cd $RESULTSDIR

# Check the value of the GENOME variable and generate individual samplesheets
if [[ "$GENOME" == "hg19" ]]; then
    python $WF_HOME/bin/split_samplesheet.py $SAMPLESHEET $RESULTSDIR
elif [[ "$GENOME" == "mm39" ]]; then
    cp $SAMPLESHEET $RESULTSDIR/mouse_rnaseq.csv
else
    echo "Unknown genome: $GENOME"
    exit 1
fi

nf_cmd="nextflow"
nf_cmd="$nf_cmd run"
nf_cmd="$nf_cmd -c $CONFIG_FILE"
nf_cmd="$nf_cmd -profile $PROFILE"
nf_cmd="$nf_cmd --logdir $LOG"
nf_cmd="$nf_cmd $WF_HOME/main.nf -resume "
nf_cmd="$nf_cmd -with-trace"
nf_cmd="$nf_cmd -with-timeline"
nf_cmd="$nf_cmd -with-dag"

echo $nf_cmd

eval $nf_cmd
