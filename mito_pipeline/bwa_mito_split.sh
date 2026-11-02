#!/bin/bash
# Usage: bwa_mito_split.sh <Out> <Job> <sample>
# Aligns to the rotated chrM so reads spanning the 16569/1 junction map properly.
# Only the +-200 bp around the junction are used from this alignment.
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

module load SAMtools/1.18-GCC-12.3.0
module load BWA/0.7.18-GCCcore-13.3.0
set -eo pipefail

read -r _ r1 r2 < <(awk -v s="$name" '$1==s' "${Out}/samplesheet.txt")
[[ -n "$r1" && -n "$r2" ]] || { echo "ERROR: ${name} not found in samplesheet" >&2; exit 1; }

raw=${Out}/${Job}_rawfiles
bwa_folder=${Out}/${Job}_results/${Job}_bam

bwa mem -t "${THREADS}" -T 19 -B 3 -L 5,4 \
    -R "@RG\tID:${name}\tSM:${name}\tPL:ILLUMINA" \
    "${SPLIT_REF}" "${raw}/${r1}" "${raw}/${r2}" \
    | samtools view -b -F 4 - \
    | samtools sort -@ "${THREADS}" -o "${bwa_folder}/${name}_split_bwa.bam" -

samtools index "${bwa_folder}/${name}_split_bwa.bam"

samtools flagstat "${bwa_folder}/${name}_split_bwa.bam" > "${bwa_folder}/${name}_stat_split_bwa.txt"
samtools idxstats "${bwa_folder}/${name}_split_bwa.bam" > "${bwa_folder}/${name}_idxstat_split_bwa.txt"

samtools view -b -q 1 "${bwa_folder}/${name}_split_bwa.bam" > "${bwa_folder}/${name}_split_bwa_accepted.bam"
samtools index "${bwa_folder}/${name}_split_bwa_accepted.bam"

echo "Split-reference alignment done for ${name}"
