#!/bin/bash
# Usage: bash folder_setup_mito.sh <input_dir> <output_dir> <job_name>
#
# Handles paired-end files named like:
#   SAMPLE_L006_R1_001.fastq(.gz)   (lane-split; lanes are concatenated)
#   SAMPLE_R1.fastq(.gz)            (already one file per read)
# Output: ${Job}_rawfiles/SAMPLE_R1.fastq.gz and SAMPLE_R2.fastq.gz
# Samplesheet: <sample> <R1 file> <R2 file>

set -euo pipefail
shopt -s nullglob

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <job_name>" >&2
    exit 1
fi

In=$(realpath "$1")
mkdir -p "$2"
Out=$(realpath "$2")
Job=$3

# Use pigz (parallel gzip) if available, otherwise gzip
if command -v pigz >/dev/null 2>&1; then GZ="pigz -c"; else GZ="gzip -c"; fi

# ---- Folders ---------------------------------------------------------------
cd "$Out"
mkdir -p "${Job}_logs" "${Job}_results" "${Job}_rawfiles"

raw=${Out}/${Job}_rawfiles
fastq_folder=${Out}/${Job}_results/${Job}_fastqc
bwa_folder=${Out}/${Job}_results/${Job}_bam
vcf_folder=${Out}/${Job}_results/${Job}_vcf
cov_folder=${Out}/${Job}_results/${Job}_cov
snp_folder=${Out}/${Job}_results/${Job}_snpEff
mkdir -p "$fastq_folder" "$bwa_folder" "$vcf_folder" "$cov_folder" "$snp_folder"

echo "Setup of the folders is done"

# ---- Find samples ----------------------------------------------------------
cd "$In"

r1_files=( *_R1_001.fastq *_R1_001.fastq.gz *_R1.fastq *_R1.fastq.gz )
if (( ${#r1_files[@]} == 0 )); then
    echo "ERROR: no R1 fastq files found in $In" >&2
    exit 1
fi

# Strip lane / read / suffix to get the sample name, e.g.
# GH_1_GTATCGTC_L006_R1_001.fastq -> GH_1_GTATCGTC
mapfile -t samples < <(printf '%s\n' "${r1_files[@]}" \
    | sed -E 's/(_L[0-9]+)?_R1(_001)?\.fastq(\.gz)?$//' | sort -u)

echo "Found ${#samples[@]} sample(s): ${samples[*]}"

# ---- Combine lanes / compress into R1 and R2 per sample --------------------
for sample in "${samples[@]}"; do
    for read in R1 R2; do
        # Lane-split files first (sorted by lane so R1 and R2 orders match),
        # otherwise single per-read files
        parts=( "${sample}"_L[0-9]*_${read}_001.fastq.gz "${sample}"_L[0-9]*_${read}_001.fastq )
        if (( ${#parts[@]} == 0 )); then
            parts=( "${sample}_${read}.fastq.gz" "${sample}_${read}.fastq" \
                    "${sample}_${read}_001.fastq.gz" "${sample}_${read}_001.fastq" )
        fi

        if (( ${#parts[@]} == 0 )); then
            echo "ERROR: no ${read} file found for sample ${sample}" >&2
            exit 1
        fi

        out=${raw}/${sample}_${read}.fastq.gz
        if [[ -s "$out" ]]; then
            echo "  ${sample} ${read}: already exists, skipping"
            continue
        fi

        echo "  ${sample} ${read}: ${parts[*]}"
        # Concatenated gzip files are valid gzip, so .gz parts can be cat'ed;
        # uncompressed parts are compressed on the fly
        {
            for p in "${parts[@]}"; do
                if [[ "$p" == *.gz ]]; then cat "$p"; else $GZ "$p"; fi
            done
        } > "${out}.tmp"
        mv "${out}.tmp" "$out"
    done
done

echo "The fastq files are combined and copied to the working directory"

# ---- Samplesheets ----------------------------------------------------------
cd "$raw"

ls -d "$PWD"/*.fastq.gz > "${Out}/sample_list.txt"

: > "${Out}/samplesheet.txt"
for sample in "${samples[@]}"; do
    echo "${sample} ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz" >> "${Out}/samplesheet.txt"
done

echo "Sample sheet generation is done:"
cat "${Out}/samplesheet.txt"