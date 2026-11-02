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

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch master_mtDNA.sh <fastq_dir> <output_dir> <job_name>" >&2
    exit 1
fi
In=$1
Out=$2
Job=$3

# ---- Find the pipeline folder ----------------------------------------------
# Under sbatch this script runs from a copy in /tmp, so ask SLURM where the
# original is; fall back to the folder sbatch was run from.
self=""
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    self=$(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null | awk -F'Command=' '/Command=/{print $2; exit}' | awk '{print $1}')
fi
[[ -z "${self}" || ! -f "${self}" ]] && self=$0
script_path=$(dirname "$(realpath "${self}")")
if [[ ! -f "${script_path}/config_mito.sh" ]]; then
    script_path=$(realpath "${SLURM_SUBMIT_DIR:-.}")
fi
if [[ ! -f "${script_path}/config_mito.sh" ]]; then
    echo "ERROR: cannot find the pipeline folder (config_mito.sh)." >&2
    echo "Submit from inside the pipeline folder, or with its full path." >&2
    exit 1
fi
source "${script_path}/config_mito.sh"

echo "Starting pipeline: ${Job}"
echo "Pipeline folder:   ${PIPE_DIR}"

# ---- Check the pipeline folder has been prepared -----------------------------
if ! check_resources; then
    echo "ERROR: the pipeline folder is not prepared (files listed above are missing)." >&2
    echo "On a login node run:  bash ${PIPE_DIR}/prepare_mito_pipeline.sh" >&2
    exit 1
fi
echo "Resources OK"

# ---- Step 1: Setup -----------------------------------------------------------
echo "Step 1: Setup"
mkdir -p "$Out"
bash "${PIPE_DIR}/folder_setup_mito.sh" "$In" "$Out" "$Job" || { echo "ERROR: setup failed"; exit 1; }
Out=$(realpath "$Out")
logs=${Out}/${Job}_logs

no_of_samples=$(wc -l < "$Out/samplesheet.txt")
echo "Found $no_of_samples samples"

# Split the job's CPUs between samples running at the same time
export THREADS=$(( ${SLURM_CPUS_PER_TASK:-1} / no_of_samples ))
(( THREADS < 1 )) && THREADS=1
echo "Using $THREADS threads per sample"

# Run one step script for every sample in parallel, one log file per sample,
# and stop the pipeline if any sample fails.
run_step () {
    local step_name=$1 script=$2
    local pids=() samples=() fail=0
    echo "${step_name}"
    while read -r sample _; do
        bash "${PIPE_DIR}/${script}" "$Out" "$Job" "$sample" \
            > "${logs}/${Job}_${sample}_${script%.sh}.log" 2>&1 < /dev/null &
        pids+=($!); samples+=("$sample")
    done < "$Out/samplesheet.txt"

    for i in "${!pids[@]}"; do
        if ! wait "${pids[$i]}"; then
            echo "  FAILED: ${samples[$i]} (see ${logs}/${Job}_${samples[$i]}_${script%.sh}.log)"
            fail=1
        fi
    done
    (( fail )) && { echo "ERROR: ${step_name} failed, stopping"; exit 1; }
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
