#!/bin/bash

set -e

# load singularity and nextflow modules
module load singularity nextflow


# set workDir ... by default it goes to `pwd`/work
# this can also be set using "workDir" in nextflow.config
# export NXF_WORK="/data/khanlab3/kopardevn/AWS_MVP_test/work"
export OUTDIR="/data/khanlab3/kopardevn/AWS_MVP_test"
export OUTTAG="3" # workDir will be $OUTDIR/work.$OUTTAG and resultsDir will be $OUTDIR/results.$OUTTAG and singularity cache is set to $OUTDIR/.singularity

# set .nextflow dir ... dont want this to go to $HOME/.nextflow
export NXF_HOME="$OUTDIR/results.$OUTTAG/.nextflow"

printenv|grep NXF

# run
nextflow run -profile biowulf main.nf -resume
