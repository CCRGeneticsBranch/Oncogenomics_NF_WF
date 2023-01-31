#!/bin/bash

if [[ "$#" -ne "1" ]];then
	echo "TAG argument is required!"
	exit
fi

# list of profiles
# biowulf_test_run_local -> get interactive node and run there
# biowulf_test_run_slurm -> get interactive node and submit jobs to slurm
# 

# PROFILE="biowulf_test_run_local"
PROFILE="biowulf_test_run_slurm"

set -e

SCRIPT_NAME="$0"
SCRIPT_DIRNAME=$(readlink -f $(dirname $0))
SCRIPT_BASENAME=$(basename $0)
WF_HOME=$SCRIPT_DIRNAME

CONFIG_FILE="$WF_HOME/nextflow.config"

# load singularity and nextflow modules
module load singularity nextflow graphviz

# set workDir ... by default it goes to `pwd`/work
# this can also be set using "workDir" in nextflow.config
# export NXF_WORK="/data/khanlab2/kopardevn/AWS_MVP_test/work"
export OUTDIR="/data/khanlab3/kopardevn/AWS_MVP_test"
# export OUTTAG="9" # workDir will be $OUTDIR/work.$OUTTAG and resultsDir will be $OUTDIR/results.$OUTTAG and singularity cache is set to $OUTDIR/.singularity
export OUTTAG=$1
export RESULTSDIR="$OUTDIR/results.$OUTTAG"
export WORKDIR="$OUTDIR/work.$OUTTAG"

# set .nextflow dir ... dont want this to go to $HOME/.nextflow
export NXF_HOME="$RESULTSDIR/.nextflow"
# export SINGULARITY_BIND="/lscratch/$SLURM_JOB_ID"

printenv|grep NXF

# run
if [ ! -d $RESULTSDIR ]; then mkdir -p $RESULTSDIR;fi
cd $RESULTSDIR
#nextflow run -profile biowulf main.nf -resume
nf_cmd="nextflow"
nf_cmd="$nf_cmd run"
nf_cmd="$nf_cmd -c $CONFIG_FILE"
nf_cmd="$nf_cmd -profile $PROFILE"
nf_cmd="$nf_cmd $WF_HOME/main.nf -resume "
# nf_cmd="$nf_cmd -with-report $RESULTSDIR/report.html"
nf_cmd="$nf_cmd -with-trace"
nf_cmd="$nf_cmd -with-timeline $RESULTSDIR/timeline.html"
# nf_cmd="$nf_cmd -with-dag $RESULTSDIR/dag.png"

echo $nf_cmd

eval $nf_cmd
