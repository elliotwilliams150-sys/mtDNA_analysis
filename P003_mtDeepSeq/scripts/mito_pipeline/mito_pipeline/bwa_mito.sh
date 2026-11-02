#!/bin/bash
# Usage: bwa_mito.sh <Out> <Job> <sample>
# Paired-end alignment to the whole GRCh38 genome (so NUMT reads go to their
# nuclear location), then keeps uniquely mapped reads on chrM.
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

# module load SAMtools/1.18-GCC-12.3.0
# module load BWA/0.7.18-GCCcore-13.3.0
set -eo pipefail

read -r _ r1 r2 < <(awk -v s="$name" '$1==s' "${Out}/samplesheet.txt")
[[ -n "$r1" && -n "$r2" ]] || { echo "ERROR: ${name} not found in samplesheet" >&2; exit 1; }

raw=${Out}/${Job}_rawfiles
bwa_folder=${Out}/${Job}_results/${Job}_bam

# -P removed: it disables read pairing, which is wanted for paired-end data
bwa mem -t "${THREADS}" -T 19 -B 3 -L 5,4 \
    -R "@RG\tID:${name}\tSM:${name}\tPL:ILLUMINA" \
    "${REF}" "${raw}/${r1}" "${raw}/${r2}" \
    | samtools sort -@ "${THREADS}" -o "${bwa_folder}/${name}_bwa.bam" -

samtools index "${bwa_folder}/${name}_bwa.bam"

# QC statistics (idxstats shows how many reads hit chrM vs nuclear contigs)
samtools flagstat "${bwa_folder}/${name}_bwa.bam" > "${bwa_folder}/${name}_stat_bwa.txt"
samtools idxstats "${bwa_folder}/${name}_bwa.bam" > "${bwa_folder}/${name}_idxstat_bwa.txt"

echo "The reads have been aligned"

# Keep uniquely mapped reads (MAPQ >= 1) on chrM only
samtools view -b -q 1 "${bwa_folder}/${name}_bwa.bam" chrM > "${bwa_folder}/${name}_bwa_accepted.bam"
samtools index "${bwa_folder}/${name}_bwa_accepted.bam"

echo "chrM reads extracted for ${name}"
