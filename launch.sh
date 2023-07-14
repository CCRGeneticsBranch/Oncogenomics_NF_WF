#!/bin/bash

if [[ "$#" -ne "2" ]];then
	echo "This script takes two inputs"
	echo "Provide Tag argument - This will add a tag to your resultsdir."
	echo "Provide samplesheet with complete path"
	exit
fi

export TAG=$1
export SAMPLESHEET=$2

sh /data/khanlab/projects/Nextflow_dev/AWS_POC_Nextflow/nf.sh $TAG  $SAMPLESHEET
