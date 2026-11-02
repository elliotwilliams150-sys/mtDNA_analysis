#!/bin/bash
# Shared settings for the human mtDNA pipeline. Sourced by every script.
# The pipeline folder is wherever this file is, so the folder can be moved
# or copied without editing any paths.

PIPE_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

# ---- Reference genome (downloaded by prepare_mito_pipeline.sh) -------------
REF_DIR=${PIPE_DIR}/ReferenceGenome
REF_NAME=GCA_000001405.15_GRCh38_no_alt_analysis_set
REF=${REF_DIR}/${REF_NAME}.fa
REF_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

# Rotated chrM reference for the 16569/1 junction
SPLIT_DIR=${REF_DIR}/split_reference
SPLIT_REF=${SPLIT_DIR}/reference_split.fa
SPLIT_DATA=${SPLIT_DIR}/reference_split.data
CHRM_BED=${SPLIT_DIR}/chrM.bed

# ---- Programs (installed by prepare_mito_pipeline.sh) ----------------------
PROGRAMS_DIR=${PIPE_DIR}/Programs
LOFREQ_SIF=${PROGRAMS_DIR}/lofreq.sif
LOFREQ_IMAGE=docker://quay.io/biocontainers/lofreq:2.1.5--py310h8360dc1_7
SNPEFF_DIR=${PROGRAMS_DIR}/snpEff
SNPEFF_URL=https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip

# Newer snpEff zips put SnpSift in its own folder; older ones inside snpEff/
if [[ -s ${PROGRAMS_DIR}/snpSift/SnpSift.jar ]]; then
    SNPSIFT_JAR=${PROGRAMS_DIR}/snpSift/SnpSift.jar
else
    SNPSIFT_JAR=${SNPEFF_DIR}/SnpSift.jar
fi

# snpEff database: Ensembl GRCh38 (mitochondrion is called "MT",
# genes labelled Mt_tRNA / Mt_rRNA / protein_coding)
SNPEFF_DB=GRCh38.115
SNPEFF_CHROM=MT

# ---- MitoSAlt: large deletions / duplications (step 8) ---------------------
RUN_MITOSALT=yes                 # set to "no" to skip step 8
MSALT_DIR=${PROGRAMS_DIR}/MitoSAlt
MSALT_URL=https://github.com/larsson-lab/MitoSAlt/archive/refs/heads/master.zip
MSALT_SCRIPT=${MSALT_DIR}/MitoSAlt1.1.1.pl
MSALT_CONFIG=${MSALT_DIR}/config_human_GRCh38.txt
MSALT_RLIB=${MSALT_DIR}/Rlib
MSALT_LAST_URL=https://gitlab.com/mcfrith/last/-/archive/1650/last-1650.tar.gz
MSALT_BBMAP_URL=https://downloads.sourceforge.net/project/bbmap/BBMap_38.05.tar.gz
MSALT_BG2BW_URL=https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
R_MODULE=R/4.4.2-gfbf-2024a     # module that provides R (check with: module spider R)

# ---- Threads per sample (master_mtDNA.sh exports THREADS) -----------------
THREADS=${THREADS:-4}

# ---- Check that everything prepare_mito_pipeline.sh creates is present ----
check_resources () {
    local missing=0 f
    for f in "${REF}" "${REF}.fai" "${REF}.amb" "${REF}.ann" "${REF}.bwt" "${REF}.pac" "${REF}.sa" \
             "${SPLIT_REF}" "${SPLIT_REF}.fai" "${SPLIT_REF}.bwt" "${SPLIT_DATA}" "${CHRM_BED}" \
             "${LOFREQ_SIF}" "${SNPEFF_DIR}/snpEff.jar" "${SNPEFF_DIR}/snpEff.config" "${SNPSIFT_JAR}" \
             "${SNPEFF_DIR}/data/${SNPEFF_DB}/snpEffectPredictor.bin"; do
        if [[ ! -s "$f" ]]; then echo "  MISSING: $f" >&2; missing=1; fi
    done
    if [[ "${RUN_MITOSALT}" == "yes" ]]; then
        for f in "${MSALT_SCRIPT}" "${MSALT_CONFIG}" "${MSALT_DIR}/delplot.R" \
                 "${MSALT_DIR}/bin/last/bin/lastal" "${MSALT_DIR}/bin/last/bin/last-split" \
                 "${MSALT_DIR}/bin/last/bin/maf-convert" "${MSALT_DIR}/bin/bbmap/reformat.sh" \
                 "${MSALT_DIR}/bin/bedGraphToBigWig" "${MSALT_DIR}/genome/human_mt_rCRS.fasta.fai" \
                 "${MSALT_DIR}/genome/human_mt_rCRS.prj" "${MSALT_DIR}/wrappers/R" \
                 "${MSALT_RLIB}/Biostrings/DESCRIPTION" "${MSALT_RLIB}/pwalign/DESCRIPTION"; do
            if [[ ! -s "$f" ]]; then echo "  MISSING (MitoSAlt): $f" >&2; missing=1; fi
        done
    fi
    return $missing
}
