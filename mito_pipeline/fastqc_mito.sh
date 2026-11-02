#!/bin/bash
# Usage: fastqc_mito.sh <Out> <Job> <sample>
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

module load FastQC/0.11.9-Java-11
set -eo pipefail

read -r _ r1 r2 < <(awk -v s="$name" '$1==s' "${Out}/samplesheet.txt")
raw=${Out}/${Job}_rawfiles
fastq_folder=${Out}/${Job}_results/${Job}_fastqc

fastqc -t 2 --outdir="${fastq_folder}" "${raw}/${r1}" "${raw}/${r2}"

echo "FastQC done for ${name}"
