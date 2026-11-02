#!/bin/bash
# Usage: lofreq_mito.sh <Out> <Job> <sample>
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

module load SAMtools/1.18-GCC-12.3.0
module load apptainer
set -eo pipefail
trap 'echo "ERROR in lofreq_mito.sh line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# Stop with a clear message if a step did not produce its output file
need () { [[ -s "$1" ]] || { echo "ERROR: expected output not created: $1" >&2; exit 1; }; }

bwa_folder=${Out}/${Job}_results/${Job}_bam
vcf_folder=${Out}/${Job}_results/${Job}_vcf
cd "${vcf_folder}"

# Make the pipeline folder and the output folder visible inside the container
BINDS="${PIPE_DIR}"
case "$(realpath "${Out}")/" in "${PIPE_DIR}"/*) ;; *) BINDS="${BINDS},$(realpath "${Out}")" ;; esac
LOFREQ="apptainer exec --bind ${BINDS} ${LOFREQ_SIF} lofreq"

# lofreq call-parallel needs at least 2 threads
lofreq_call () {
    if (( THREADS > 1 )); then
        $LOFREQ call-parallel --pp-threads "${THREADS}" "$@"
    else
        $LOFREQ call "$@"
    fi
}

# --------------------------------------------------
# Indel qualities (intermediate BAMs, kept in the BAM folder)
# --------------------------------------------------
iq_main=${bwa_folder}/${name}_indelqual.bam
iq_junc=${bwa_folder}/${name}_junction_indelqual.bam
rm -f ${iq_main}* ${iq_junc}* ${name}_indelqual.bam* ${name}_junction_indelqual.bam*

$LOFREQ indelqual --dindel --ref "${REF}" \
    --out ${iq_main} ${bwa_folder}/${name}_bwa_accepted.bam
$LOFREQ indelqual --dindel --ref "${SPLIT_REF}" \
    --out ${iq_junc} ${bwa_folder}/${name}_split_bwa_accepted.bam
need ${iq_main}; need ${iq_junc}

samtools index ${iq_main}
samtools index ${iq_junc}

# --------------------------------------------------
# Variant calling (main alignment restricted to chrM)
# --------------------------------------------------
rm -f ${name}_nofilter.vcf ${name}_junction_nofilter.vcf

lofreq_call -f "${REF}" -l "${CHRM_BED}" \
    -o ${name}_nofilter.vcf \
    -N -B -q 30 -Q 30 --call-indels --no-default-filter \
    ${iq_main}
need ${name}_nofilter.vcf

lofreq_call -f "${SPLIT_REF}" \
    -o ${name}_junction_nofilter.vcf \
    -N -B -q 30 -Q 30 --call-indels --no-default-filter \
    ${iq_junc}
need ${name}_junction_nofilter.vcf

# --------------------------------------------------
# Filtering
# --------------------------------------------------
rm -f ${name}_sbfiltered.vcf ${name}_junction_sbfiltered.vcf

$LOFREQ filter --no-defaults --snvqual-thresh 70 --indelqual-thresh 70 \
    --sb-incl-indels -B 60 -i ${name}_nofilter.vcf -o ${name}_sbfiltered.vcf

$LOFREQ filter --no-defaults --snvqual-thresh 70 --indelqual-thresh 70 \
    --sb-incl-indels -B 60 -i ${name}_junction_nofilter.vcf -o ${name}_junction_sbfiltered.vcf
need ${name}_sbfiltered.vcf; need ${name}_junction_sbfiltered.vcf

# --------------------------------------------------
# Junction processing: replace +-200 bp around 16569/1 with split-reference calls
# --------------------------------------------------
mid_point=$(awk '{print $2}' "${SPLIT_DATA}")
end_point=$(awk '{print $3}' "${SPLIT_DATA}")

awk '/#CHROM/ {flag = 1; next} flag {print}' ${name}_junction_sbfiltered.vcf | \
awk -v midpoint="$mid_point" -v endpoint="$end_point" \
'BEGIN {OFS = "\t"}
{test = $2 - midpoint; if (test <= 0) res = endpoint + test; else res = test; $2 = res; print}' | \
sort -nk2 | \
awk -v endpoint="$end_point" '$2 <= 200 || $2 >= endpoint - 200 {print}' \
> ${name}_junction_replacement.vars

awk '/^#/ {print}' ${name}_sbfiltered.vcf > ${name}_sbfiltered_junction_combined.vcf

awk '/#CHROM/ {flag = 1; next} flag {print}' ${name}_sbfiltered.vcf | \
awk -v endpoint="$end_point" '$2 > 200 && $2 < endpoint - 200 {print}' \
> ${name}_middle.vars

cat ${name}_middle.vars ${name}_junction_replacement.vars | \
sort -nk2 >> ${name}_sbfiltered_junction_combined.vcf

# --------------------------------------------------
# Split SNVs and indels
# --------------------------------------------------
rm -f ${name}_snvs.vcf ${name}_indels.vcf

$LOFREQ filter --no-defaults --only-snvs \
    -i ${name}_sbfiltered_junction_combined.vcf -o ${name}_snvs.vcf

$LOFREQ filter --no-defaults --only-indels \
    -i ${name}_sbfiltered_junction_combined.vcf -o ${name}_indels.vcf

need ${name}_snvs.vcf; need ${name}_indels.vcf

echo "Variant calling done for ${name}:"
for f in ${name}_nofilter.vcf ${name}_sbfiltered_junction_combined.vcf ${name}_snvs.vcf ${name}_indels.vcf; do
    echo "  $(grep -vc '^#' "$f" || true) variants in ${vcf_folder}/$f"
done
