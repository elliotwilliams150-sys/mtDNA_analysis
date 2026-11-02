#!/bin/bash
# Usage: mitosalt_mito.sh <Out> <Job> <sample>
# Step 8: detects large mtDNA deletions and duplications with MitoSAlt.
#
# Reads that map to chrM (plus their mates, and pairs where neither read maps)
# are taken from the whole-genome BWA alignment of step 3. These are the same
# reads MitoSAlt would select itself, so MitoSAlt is run in "enriched" mode and
# goes straight to split-read alignment against chrM with LAST.
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

module load SAMtools/1.18-GCC-12.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load Java/21.0.5
set -eo pipefail
trap 'echo "ERROR in mitosalt_mito.sh line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

need () { [[ -s "$1" ]] || { echo "ERROR: expected output not created: $1" >&2; exit 1; }; }

bwa_bam=${Out}/${Job}_results/${Job}_bam/${name}_bwa.bam
msalt_folder=${Out}/${Job}_results/${Job}_mitosalt
work=${msalt_folder}/${name}
need "${bwa_bam}"

# MitoSAlt writes into fixed sub-folders of the folder it is run from
rm -rf "${work}"
mkdir -p "${work}"/{bam,bw,indel,log,plot,tab}
cd "${work}"
ln -sf "${MSALT_DIR}/delplot.R" delplot.R

# R (for the plots) comes from its own module via the wrapper
export PATH="${MSALT_DIR}/wrappers:${PATH}"

# --------------------------------------------------
# 1. Extract mitochondrial read pairs from the BWA alignment
# --------------------------------------------------
echo "Extracting chrM read pairs for ${name}"
samtools view -@ "${THREADS}" -u -F 0x900 \
    -e 'rname=="chrM" || mrname=="chrM" || (flag.unmap && flag.munmap)' "${bwa_bam}" \
  | samtools collate -@ "${THREADS}" -u -O - "tmp_collate_${name}" \
  | samtools fastq -@ "${THREADS}" -n \
      -1 "${name}_mt_R1.fq.gz" -2 "${name}_mt_R2.fq.gz" -0 /dev/null -s /dev/null -
need "${name}_mt_R1.fq.gz"; need "${name}_mt_R2.fq.gz"
pairs=$(( $(zcat "${name}_mt_R1.fq.gz" | wc -l) / 4 ))
echo "  ${pairs} read pairs"
(( pairs > 0 )) || { echo "ERROR: no chrM reads found for ${name}" >&2; exit 1; }

# --------------------------------------------------
# 2. Run MitoSAlt
# --------------------------------------------------
sed "s/^threads = .*/threads = ${THREADS}/" "${MSALT_CONFIG}" > config_${name}.txt

perl "${MSALT_SCRIPT}" "config_${name}.txt" \
    "${name}_mt_R1.fq.gz" "${name}_mt_R2.fq.gz" "${name}"

# MitoSAlt does not stop on errors, so check its outputs here
need "bam/${name}.bam"
need "tab/${name}.tab.gz"
[[ -e "indel/${name}.cluster" ]] || { echo "ERROR: MitoSAlt did not create indel/${name}.cluster" >&2; exit 1; }
if [[ -f delplot.Rout ]] && grep -q "Execution halted" delplot.Rout; then
    tail -20 delplot.Rout >&2
    echo "ERROR: MitoSAlt's R plotting script failed (see ${work}/delplot.Rout)" >&2
    exit 1
fi

# --------------------------------------------------
# 3. Collect the main results one level up
# --------------------------------------------------
if [[ -s "indel/${name}.tsv" ]]; then
    cp "indel/${name}.tsv" "${msalt_folder}/${name}_mitosalt.tsv"
    events=$(( $(wc -l < "indel/${name}.tsv") - 1 ))
else
    # No deletion/duplication passed MitoSAlt's filters
    printf "sample\tcluster.id\talt.reads\tref.reads\theteroplasmy\tdel.start.range\tdel.end.range\tdel.size\tfinal.event\tfinal.start\tfinal.end\tfinal.size\tseq1\tseq2\tseq\n" \
        > "${msalt_folder}/${name}_mitosalt.tsv"
    events=0
fi
[[ -s "plot/${name}.pdf" ]] && cp "plot/${name}.pdf" "${msalt_folder}/${name}_mitosalt.pdf"

rm -f "${name}_mt_R1.fq.gz" "${name}_mt_R2.fq.gz"

echo "MitoSAlt done for ${name}: ${events} deletion/duplication event(s)"
echo "  table: ${msalt_folder}/${name}_mitosalt.tsv"
