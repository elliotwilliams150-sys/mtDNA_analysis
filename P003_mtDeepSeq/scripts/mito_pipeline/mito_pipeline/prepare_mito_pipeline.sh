#!/bin/bash
# Prepares the pipeline folder: downloads the reference genome and builds the
# split chrM reference (ReferenceGenome/), and installs LoFreq, snpEff/SnpSift,
# the snpEff database and MitoSAlt with its tools (Programs/).
#
# Run ONCE, on a LOGIN node (needs internet), from anywhere:
#   bash /path/to/mito_pipeline/prepare_mito_pipeline.sh
#
# Safe to re-run: anything already in place is skipped.
set -eo pipefail

source "$(dirname "$(realpath "$0")")/config_mito.sh"

#module load SAMtools/1.18-GCC-12.3.0
#module load BWA/0.7.18-GCCcore-13.3.0
#module load apptainer
#module load Java/21.0.5

# Use software available in the current server/Conda environment.
required=(samtools bwa java apptainer wget unzip tar gzip md5sum)

if [[ "${RUN_MITOSALT}" == "yes" ]]; then
    required+=(bedtools R Rscript python3 perl make g++)
fi

missing=()
for cmd in "${required[@]}"; do
    command -v "$cmd" >/dev/null 2>&1 || missing+=("$cmd")
done

if ((${#missing[@]})); then
    printf 'ERROR: missing commands: %s\n' "${missing[*]}" >&2
    exit 1
fi

java_version=$(java -version 2>&1)
java_major=$(sed -n 's/.*version "\([0-9]*\).*/\1/p' <<< "$java_version" | head -1)

if [[ ! "$java_major" =~ ^[0-9]+$ ]] || ((java_major < 21)); then
    echo "ERROR: Java 21 or newer is required for the current snpEff download." >&2
    printf '%s\n' "$java_version" >&2
    exit 1
fi


log () { echo "[$(date '+%H:%M:%S')] $*"; }

log "Preparing pipeline folder: ${PIPE_DIR}"
mkdir -p "${REF_DIR}" "${SPLIT_DIR}" "${PROGRAMS_DIR}"

# =============================================================================
# 1. Reference genome: GRCh38 no-ALT analysis set, with BWA index
# =============================================================================
log "1/6 Reference genome"
cd "${REF_DIR}"

ref_ready () {
    local e
    for e in "" .fai .amb .ann .bwt .pac .sa; do [[ -s "${REF}${e}" ]] || return 1; done
}

if ref_ready; then
    log "    already present, skipping"
else
    wget -q -c "${REF_URL}/md5checksums.txt"

    # Download one file and check its MD5
    fetch () {
        local f=$1 expected actual
        log "    downloading ${f}"
        wget -c -q --show-progress "${REF_URL}/${f}"
        expected=$(grep -E "[ /]${f}$" md5checksums.txt | awk '{print $1}' || true)
        actual=$(md5sum "${f}" | awk '{print $1}')
        if [[ -n "${expected}" && "${expected}" != "${actual}" ]]; then
            echo "ERROR: checksum mismatch for ${f}. Delete it and re-run." >&2
            exit 1
        fi
    }

    if [[ ! -s "${REF}" ]]; then
        fetch "${REF_NAME}.fna.gz"
        log "    decompressing"
        gunzip -c "${REF_NAME}.fna.gz" > "${REF}.tmp" && mv "${REF}.tmp" "${REF}"
        rm -f "${REF_NAME}.fna.gz"
    fi

    if [[ ! -s "${REF}.fai" ]]; then
        fetch "${REF_NAME}.fna.fai"
        mv "${REF_NAME}.fna.fai" "${REF}.fai"
    fi

    if ! ref_ready; then
        fetch "${REF_NAME}.fna.bwa_index.tar.gz"
        log "    extracting BWA index"
        tar -xzf "${REF_NAME}.fna.bwa_index.tar.gz"
        for f in "${REF_NAME}".fna.*; do
            [[ "$f" == *.tar.gz ]] && continue
            mv -f "$f" "${f/${REF_NAME}.fna./${REF_NAME}.fa.}"
        done
        rm -f "${REF_NAME}.fna.bwa_index.tar.gz"
    fi

    rm -f md5checksums.txt
    ref_ready || { echo "ERROR: reference files incomplete in ${REF_DIR}" >&2; exit 1; }
    log "    done"
fi

# =============================================================================
# 2. Split (rotated) chrM reference for the 16569/1 junction
# =============================================================================
log "2/6 Split chrM reference"
cd "${SPLIT_DIR}"

samtools faidx "${REF}" chrM > chrM.fa
samtools faidx chrM.fa
L=$(cut -f2 chrM.fa.fai)          # 16569 for rCRS
MID=$(( L / 2 ))

# Last MID bases followed by the first L-MID bases, so the junction sits in
# the middle. The contig stays "chrM" so junction calls merge cleanly.
seq=$(grep -v '>' chrM.fa | tr -d '\n')
{ echo ">chrM"; echo "${seq:L-MID}${seq:0:L-MID}" | fold -w 60; } > reference_split.fa

# <contig> <mid_point> <end_point>
printf "chrM\t%s\t%s\n" "${MID}" "${L}" > reference_split.data
printf "chrM\t0\t%s\n" "${L}" > chrM.bed

samtools faidx reference_split.fa
bwa index reference_split.fa 2> bwa_index.log
log "    done (chrM length ${L}, mid point ${MID})"

# =============================================================================
# 3. LoFreq container
# =============================================================================
log "3/6 LoFreq"
cd "${PROGRAMS_DIR}"
if [[ -s "${LOFREQ_SIF}" ]]; then
    log "    already present, skipping"
else
    apptainer pull "${LOFREQ_SIF}" "${LOFREQ_IMAGE}"
fi
log "    $(apptainer exec "${LOFREQ_SIF}" lofreq version 2>&1 | head -1)"

# =============================================================================
# 4. snpEff and SnpSift
# =============================================================================
log "4/6 snpEff / SnpSift"
cd "${PROGRAMS_DIR}"
if [[ -s "${SNPEFF_DIR}/snpEff.jar" ]]; then
    log "    already present, skipping"
else
    wget -c -q --show-progress "${SNPEFF_URL}" -O snpEff_latest_core.zip
    unzip -q -o snpEff_latest_core.zip
    rm -f snpEff_latest_core.zip
fi
# Re-read the config so SNPSIFT_JAR points at wherever the zip put it
source "${PIPE_DIR}/config_mito.sh"
[[ -s "${SNPSIFT_JAR}" ]] || { echo "ERROR: SnpSift.jar not found under ${PROGRAMS_DIR}" >&2; exit 1; }
log "    snpEff:  ${SNPEFF_DIR}/snpEff.jar"
log "    SnpSift: ${SNPSIFT_JAR}"

# =============================================================================
# 5. snpEff database
# =============================================================================
log "5/6 snpEff database ${SNPEFF_DB}"
cd "${SNPEFF_DIR}"
if [[ -s "data/${SNPEFF_DB}/snpEffectPredictor.bin" ]]; then
    log "    already present, skipping"
else
    if ! grep -q "^${SNPEFF_DB}.genome" snpEff.config; then
        echo "ERROR: ${SNPEFF_DB} is not offered by this snpEff version. Ensembl GRCh38 databases available:" >&2
        grep -E "^GRCh38\.[0-9]+\.genome" snpEff.config | awk -F'.genome' '{print "   "$1}' >&2
        echo "Set SNPEFF_DB in config_mito.sh to one of these and re-run." >&2
        exit 1
    fi
    java -Xmx8g -jar snpEff.jar download -c snpEff.config "${SNPEFF_DB}"
fi

# Translate mitochondrial genes with the vertebrate mitochondrial code
if ! grep -q "^${SNPEFF_DB}.${SNPEFF_CHROM}.codonTable" snpEff.config; then
    echo "${SNPEFF_DB}.${SNPEFF_CHROM}.codonTable : Vertebrate_Mitochondrial" >> snpEff.config
    log "    added mitochondrial codon table to snpEff.config"
fi

# =============================================================================
# 6. MitoSAlt (large deletions / duplications) and the tools it needs
# =============================================================================
if [[ "${RUN_MITOSALT}" == "yes" ]]; then
    log "6/6 MitoSAlt"
    mkdir -p "${MSALT_DIR}"
    cd "${MSALT_DIR}"
    mkdir -p bin genome wrappers "${MSALT_RLIB}"

    # 6a. MitoSAlt scripts from GitHub
    if [[ -s "${MSALT_SCRIPT}" ]]; then
        log "    MitoSAlt scripts already present"
    else
        log "    downloading MitoSAlt"
        wget -q -O mitosalt.zip "${MSALT_URL}"
        unzip -q -o mitosalt.zip
        cp -r MitoSAlt-master/. .
        rm -rf MitoSAlt-master mitosalt.zip
        # MitoSAlt asks Java for 100 GB of memory; 8 GB is plenty for mtDNA reads
        sed -i 's/-Xmx100g/-Xmx8g/g' MitoSAlt1.1.1.pl MitoSAlt_SE1.1.1.pl
    fi

    # 6b. LAST aligner (version used by MitoSAlt), compiled from source
    if [[ -x bin/last/bin/lastal && -x bin/last/bin/last-split ]]; then
        log "    LAST already present"
    else
        log "    downloading and compiling LAST (a few minutes)"
        rm -rf bin/last last-1650
        wget -q -O last.tar.gz "${MSALT_LAST_URL}"
        tar -xzf last.tar.gz && rm -f last.tar.gz
        if ! make -C last-1650/src -j 4 > last_build.log 2>&1; then
            tail -20 last_build.log >&2
            echo "ERROR: LAST failed to compile (full log: ${MSALT_DIR}/last_build.log)" >&2
            exit 1
        fi
        mv last-1650 bin/last
        # Make sure the compiled programs are in bin/last/bin
        for b in lastal lastdb last-split; do
            [[ -x bin/last/bin/$b ]] || cp bin/last/src/$b bin/last/bin/
        done
    fi
    command -v python3 >/dev/null || { echo "ERROR: python3 is needed by LAST's maf-convert" >&2; exit 1; }

    # 6c. BBMap (reformat.sh, used to interleave the read pairs)
    if [[ -s bin/bbmap/reformat.sh ]]; then
        log "    BBMap already present"
    else
        log "    downloading BBMap"
        wget -q -O bbmap.tar.gz "${MSALT_BBMAP_URL}"
        tar -xzf bbmap.tar.gz -C bin && rm -f bbmap.tar.gz
    fi

    # 6d. UCSC bedGraphToBigWig (coverage track)
    if [[ -x bin/bedGraphToBigWig ]]; then
        log "    bedGraphToBigWig already present"
    else
        wget -q -O bin/bedGraphToBigWig "${MSALT_BG2BW_URL}"
        chmod 755 bin/bedGraphToBigWig
    fi

    # 6e. Mitochondrial genome for MitoSAlt: chrM (rCRS) from our GRCh38 reference
    cp "${SPLIT_DIR}/chrM.fa" genome/human_mt_rCRS.fasta
    samtools faidx genome/human_mt_rCRS.fasta
    bin/last/bin/lastdb -uNEAR genome/human_mt_rCRS genome/human_mt_rCRS.fasta
    log "    LAST index of chrM built"

    # 6f. R wrapper: runs R from its own module in a clean environment, so
    #     loading R does not clash with the SAMtools/BEDTools toolchains
    cat > wrappers/R << 'WRAP'
#!/bin/bash
WRAPDIR=$(dirname "$(realpath "$0")")
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$WRAPDIR" | paste -sd: -)
if ! type module >/dev/null 2>&1 && [[ -n "${LMOD_PKG:-}" ]]; then
    source "${LMOD_PKG}/init/bash"
fi
module purge >/dev/null 2>&1
module load __R_MODULE__ >/dev/null 2>&1
command -v "$(basename "$0")" >/dev/null || { echo "R wrapper: could not load module __R_MODULE__" >&2; exit 1; }
export R_LIBS="__RLIB__"
exec "$(basename "$0")" "$@"
WRAP
    sed -i -e "s#__R_MODULE__#${R_MODULE}#g" -e "s#__RLIB__#${MSALT_RLIB}#g" wrappers/R
    cp wrappers/R wrappers/Rscript
        # 6f. Use the R installation available in the current environment.
    # Resolve the actual binaries before creating the wrappers.
    REAL_R=$(type -P R)
    REAL_RSCRIPT=$(type -P Rscript)

    [[ -x "${REAL_R}" && -x "${REAL_RSCRIPT}" ]] || {
        echo "ERROR: R and Rscript must be available on PATH." >&2
        exit 1
    }

    # Separate private packages by R major.minor version.
    R_VERSION=$("${REAL_RSCRIPT}" -e \
        'cat(paste(R.version$major, strsplit(R.version$minor, ".", fixed=TRUE)[[1]][1], sep="."))')
    R_PACKAGE_LIB="${MSALT_RLIB}/${R_VERSION}"
    mkdir -p "${R_PACKAGE_LIB}"

    {
        echo '#!/bin/bash'
        printf 'export R_LIBS=%q\n' "${R_PACKAGE_LIB}"
        printf 'exec %q "$@"\n' "${REAL_R}"
    } > wrappers/R

    {
        echo '#!/bin/bash'
        printf 'export R_LIBS=%q\n' "${R_PACKAGE_LIB}"
        printf 'exec %q "$@"\n' "${REAL_RSCRIPT}"
    } > wrappers/Rscript

    chmod 755 wrappers/R wrappers/Rscript

    # 6g. R packages used by MitoSAlt's plotting script
    log "    checking R packages (plotrix, RColorBrewer, Biostrings, pwalign); first install takes a while"
    wrappers/Rscript -e '
        lib  <- Sys.getenv("R_LIBS"); repo <- "https://cloud.r-project.org"
        have <- rownames(installed.packages())
        need <- setdiff(c("plotrix", "RColorBrewer", "BiocManager"), have)
        if (length(need)) install.packages(need, lib = lib, repos = repo)
        # pwalign: since Bioconductor 3.19 the alignment functions MitoSAlt uses
        # (nucleotideSubstitutionMatrix, pairwiseAlignment) live in pwalign
        bioc <- setdiff(c("Biostrings", "pwalign"), rownames(installed.packages()))
        if (length(bioc))
            BiocManager::install(bioc, lib = lib, update = FALSE, ask = FALSE)
        suppressMessages({ library(plotrix); library(RColorBrewer); library(Biostrings); library(pwalign) })
        cat("    R packages OK\n")'

    # 6h. MitoSAlt configuration for this pipeline
    #     Reads come from the BWA whole-genome alignment (step 3), so MitoSAlt
    #     runs in "enriched" mode: no HISAT2 re-alignment, straight to LAST on chrM.
    cat > "${MSALT_CONFIG}" << CFG
#TOOLS (written by prepare_mito_pipeline.sh)
hisat2 = none
lastal = ${MSALT_DIR}/bin/last/bin/lastal
lastsp = ${MSALT_DIR}/bin/last/bin/last-split
mfcv = ${MSALT_DIR}/bin/last/bin/maf-convert
reformat = ${MSALT_DIR}/bin/bbmap/reformat.sh
samtools = samtools
sambamba = none
b2fq = bedtools bamtofastq
gcov = bedtools genomecov
intersectBed = bedtools intersect
sortBed = bedtools sort
clusterBed = bedtools cluster
randomBed = bedtools random
groupBy = bedtools groupby
bg2bw = ${MSALT_DIR}/bin/bedGraphToBigWig

#DATABASES
hsindex = none
faindex = ${REF}.fai
lastindex = ${MSALT_DIR}/genome/human_mt_rCRS
mtfaindex = ${MSALT_DIR}/genome/human_mt_rCRS.fasta.fai
gsize = none
MT_fasta = ${MSALT_DIR}/genome/human_mt_rCRS.fasta

#COMPUTATION
threads = 4

#MITOCHONDRIA FEATURES (rCRS)
refchr = chrM
msize = 16569
exclude = 5
orihs = 16081
orihe = 407
orils = 5730
orile = 5763

#SCORING AND FILTERING FEATURES (MitoSAlt defaults)
score_threshold = 80
evalue_threshold = 0.00001
split_length = 15
paired_distance = 1000
deletion_threshold_min = 30
deletion_threshold_max = 30000
breakthreshold = -2
cluster_threshold = 5
breakspan = 15
sizelimit = 10000
hplimit = 0.01
flank = 15
split_distance_threshold = 5

#STEPS
dna = yes
enriched = yes
nu_mt = no
rmtmp = yes
o_mt = yes
i_del = yes
cn_mt = no
CFG
    log "    MitoSAlt ready"
else
    log "6/6 MitoSAlt skipped (RUN_MITOSALT=no in config_mito.sh)"
fi

# =============================================================================
# Final check
# =============================================================================
if check_resources; then
    log "Preparation complete."
    log "Before running samples, adapt master_mtDNA.sh and the step scripts for this server."
    
else
    echo "ERROR: some files are still missing (listed above)." >&2
    exit 1
fi
