#!/bin/bash
#SBATCH --job-name=mtDNA
#SBATCH --account=rockhpc_ghwtcmr
#SBATCH --partition=default_free
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#
# Runs the mtDNA pipeline on a folder of paired-end FASTQ files.
# Usage: sbatch /path/to/mito_pipeline/master_mtDNA.sh <fastq_dir> <output_dir> <job_name>
#
# prepare_mito_pipeline.sh must have been run once beforehand.


set -eo pipefail 

if [[ $# -ne 3 ]]; then
    echo "Usage: bash master_mtDNA.sh <fastq_dir> <output_dir> <job_name>" >&2
    exit 1
fi

In=$(realpath "$1")
Out=$(realpath -m "$2")
Job=$3

# Threads used by each step script.
export THREADS="${THREADS:-4}"

if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer." >&2
    exit 1
fi

# ---- Find the pipeline folder ----------------------------------------------
# Under sbatch this script runs from a copy in /tmp, so ask SLURM where the
# original is; fall back to the folder sbatch was run from.

script_path=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

if [[ ! -f "${script_path}/config_mito.sh" ]]; then
    echo "ERROR: config_mito.sh not found beside this script." >&2
    exit 1
fi

source "${script_path}/config_mito.sh"



# ---- Check the pipeline folder has been prepared -----------------------------
if ! check_resources; then
    echo "ERROR: the pipeline folder is not prepared (files listed above are missing)." >&2
    echo "Run: bash ${PIPE_DIR}/prepare_mito_pipeline.sh" >&2
    exit 1
fi
echo "Resources OK"

# ---- Step 1: Setup -----------------------------------------------------------
echo "Step 1: Setup"
mkdir -p "$Out"
bash "${PIPE_DIR}/folder_setup_mito.sh" "$In" "$Out" "$Job" || { echo "ERROR: setup failed"; exit 1; }
Out=$(realpath "$Out")
logs=${Out}/${Job}_logs

mkdir -p "$logs"

no_of_samples=$(awk 'NF {n++} END {print n+0}' "$Out/samplesheet.txt")

if (( no_of_samples == 0 )); then
    echo "ERROR: no samples found in $Out/samplesheet.txt" >&2
    exit 1
fi

echo "Found $no_of_samples samples"
echo "Running one sample at a time, THREADS=$THREADS"

run_step () {
    local step_name=$1 script=$2
    local sample logfile

    echo "$step_name"

    while read -r sample _; do
        [[ -n "$sample" ]] || continue

        logfile="${logs}/${Job}_${sample}_${script%.sh}.log"
        echo "  Running: $sample"

        if ! bash "${PIPE_DIR}/${script}" "$Out" "$Job" "$sample" \
            > "$logfile" 2>&1 < /dev/null; then
            echo "ERROR: ${step_name} failed for ${sample}" >&2
            echo "Log: $logfile" >&2
            tail -n 30 "$logfile" >&2
            exit 1
        fi
    done < "$Out/samplesheet.txt"

    echo "  done"
}



run_step "Step 2: FastQC"              fastqc_mito.sh
run_step "Step 3: BWA alignment"       bwa_mito.sh
run_step "Step 4: BWA split reference" bwa_mito_split.sh
run_step "Step 5: Coverage"            coverage_mtDNA.sh
run_step "Step 6: Variant calling"     lofreq_mito.sh
run_step "Step 7: Annotation"          snpeff_mito.sh
if [[ "${RUN_MITOSALT}" == "yes" ]]; then
    run_step "Step 8: MitoSAlt (large deletions/duplications)" mitosalt_mito.sh
else
    echo "Step 8: MitoSAlt skipped (RUN_MITOSALT=no in config_mito.sh)"
fi

echo "Pipeline complete: ${Job}"
echo "Final annotated tables: ${Out}/${Job}_results/${Job}_snpEff/*_ann_MT.vcf"
if [[ "${RUN_MITOSALT}" == "yes" ]]; then
    echo "Deletions/duplications: ${Out}/${Job}_results/${Job}_mitosalt/*_mitosalt.tsv"
fi

