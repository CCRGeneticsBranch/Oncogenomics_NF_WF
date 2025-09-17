#!/bin/bash

# Set default values
DEFAULT_OUTDIR="/data/khanlab/projects/processed_DATA"
DEFAULT_GENOME="hg19"
DEFAULT_PLATFORM="biowulf"

# Detect and strip a bare "cleanup" arg
DO_CLEANUP=false
POSITIONAL=()
for arg in "$@"; do
  if [[ "$arg" == "cleanup" ]]; then
    DO_CLEANUP=true
  else
    POSITIONAL+=("$arg")
  fi
done
set -- "${POSITIONAL[@]}"


# Help/usage
usage() {
  cat <<EOF
Usage:
  $(basename "$0") <samplesheet.csv> [outdir] [genome] [platform] [cleanup]

Arguments:
  samplesheet.csv   Required. Path to a CSV with headers including: sample, casename.
  outdir            Optional. Results directory. Default: ${DEFAULT_OUTDIR}
  genome            Optional. One of: hg19 | mm39. Default: ${DEFAULT_GENOME}
  platform          Optional. Free-form string passed to --platform. Default: ${DEFAULT_PLATFORM}
  cleanup           Optional. Literal word "cleanup" (any position). If present, and the job
                    finishes successfully, the script will run:
                      - nextflow clean -f "<run_name>"
                      - rm -rf .nextflow

Behavior:
  • Submits a SLURM job via sbatch, writing logs to:
      <outdir>/<PATIENT>/<CASENAME>/<logname>_<SLURM_JOB_ID>_<timestamp}.out
  • Chooses Nextflow profile by genome:
      hg19 -> biowulf_test_run_slurm
      mm39 -> biowulf_mouse_RNA_slurm
  • Passes through to Nextflow:
      --samplesheet <samplesheet.csv> --resultsdir <outdir> --genome_v <genome> --platform <platform>

Examples:
  # Minimal (defaults outdir, genome=hg19, platform=biowulf)
  $(basename "$0") /path/to/merged_samplesheet.csv

  # Specify outdir and genome
  $(basename "$0") /path/to/merged_samplesheet.csv /data/khanlab/projects/processed_DATA hg19

  # Specify platform and enable cleanup
  $(basename "$0") /path/sheet.csv /data/out hg19 biowulf cleanup

EOF
}

# Show help if asked
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

# Require at least the samplesheet
if [[ "$#" -lt 1 ]]; then
  usage
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
CONFIG_FILE="$WF_HOME/nextflow.config"

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
RUN_NAME="${logname}_${timestamp}"

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name="$logname"
#SBATCH --output="$OUTDIR/$PATIENT/$CASENAME/${logname}_%A_${timestamp}.out"
#SBATCH --cpus-per-task=2
#SBATCH --mem=05g
#SBATCH --time=08-00:00:00

logfile="$OUTDIR/$PATIENT/$CASENAME/${logname}_${SLURM_JOB_ID}_${timestamp}.out"


module load nextflow/23.10.0 singularity graphviz
nextflow run -c $CONFIG_FILE -profile $PROFILE --logdir $LOG $WF_HOME/main.nf -resume --samplesheet $SAMPLESHEET --resultsdir $OUTDIR --genome_v $GENOME --platform $PLATFORM -name "$RUN_NAME"

if [[ $DO_CLEANUP == true && -f "\$logfile" ]]; then
  if tail -n 200 "\$logfile" | grep -q "Completed at:" && \
     tail -n 200 "\$logfile" | grep -q "Succeeded"; then
    echo "✅ Run succeeded — cleaning up work dirs for $RUN_NAME"
    nextflow clean -f "$RUN_NAME"
    rm -rf .nextflow
  else
    echo "⚠️  Run did not finish — skipping cleanup"
  fi
fi

EOT
