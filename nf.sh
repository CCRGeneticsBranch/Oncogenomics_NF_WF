#!/bin/bash

set -e
module load singularity nextflow
nextflow run -profile biowulf main.nf -resume
