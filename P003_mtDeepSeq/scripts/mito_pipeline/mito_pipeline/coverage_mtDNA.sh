#!/bin/bash
# Usage: coverage_mtDNA.sh <Out> <Job> <sample>
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

# module load SAMtools/1.18-GCC-12.3.0
set -eo pipefail

bwa_folder=${Out}/${Job}_results/${Job}_bam
cov_folder=${Out}/${Job}_results/${Job}_cov
cd "${cov_folder}"

# Per-base depth on chrM only (same 3 columns as bedtools genomecov -d,
# but avoids writing a 3-billion-line file for the whole genome)
samtools depth -a -r chrM "${bwa_folder}/${name}_bwa_accepted.bam" > "${name}_coverage.txt"
samtools depth -a "${bwa_folder}/${name}_split_bwa_accepted.bam" > "${name}_coverage_junction.txt"

## Combine the coverages to represent the entire mtDNA genome
mid_point=$(awk '{print $2}' "${SPLIT_DATA}")
end_point=$(awk '{print $3}' "${SPLIT_DATA}")

# Middle part from the normal alignment
awk -v endpoint="$end_point" '$2 > 200 && $2 < endpoint - 200 {print}' \
    "${name}_coverage.txt" > "${name}_coverage.middle"

# Re-coordinate the split alignment and keep only the +-200 bp around the junction
awk -v midpoint="$mid_point" -v endpoint="$end_point" \
 'BEGIN {OFS = "\t"}
 {test = $2 - midpoint; if (test <= 0) res = endpoint + test; else res = test; $2 = res; print}' \
 "${name}_coverage_junction.txt" | sort -nk2 | \
 awk -v endpoint="$end_point" '$2 <= 200 || $2 >= endpoint - 200 {print}' \
 > "${name}_coverage_junction_replacement"

cat "${name}_coverage.middle" "${name}_coverage_junction_replacement" | sort -nk2 > "${name}_coverage_final.txt"

echo "Coverage done for ${name}"
