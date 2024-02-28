#!/bin/bash


if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 <samplesheet_with_full_path> <output_directory>"
    echo "This script requires two positional arguments:"
    echo "1. Path to samplesheet"
    echo "2. Path to results directory"
    exit 1
fi


export SAMPLESHEET=$1
export OUTDIR=$2
export PATIENT=$(basename "$SAMPLESHEET" | cut -d'_' -f1)
export CASENAME=$(basename "$SAMPLESHEET" | cut -d"_" -f2,3|sed s'/.csv//g')
logname=$(basename "$SAMPLESHEET" .csv)
timestamp=$(date +"%Y%m%d-%H%M%S")
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name="$logname"
#SBATCH --output="$logname"_%A_$timestamp.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=50g
#SBATCH --time=10:00:00


sh /data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/nf.sh $SAMPLESHEET $OUTDIR $PATIENT $CASENAME

exit 0
EOT
