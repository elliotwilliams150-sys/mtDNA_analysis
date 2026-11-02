#!/bin/bash
# Usage: snpeff_mito.sh <Out> <Job> <sample>
Out=$1; Job=$2; name=$3

source "$(dirname "$(realpath "$0")")/config_mito.sh"

module load Java/21.0.5
set -eo pipefail
trap 'echo "ERROR in snpeff_mito.sh line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

lofreq_dir=${Out}/${Job}_results/${Job}_vcf
snpeff_dir=${Out}/${Job}_results/${Job}_snpEff
cd "${snpeff_dir}"

# ---- Pre-flight checks ----------------------------------------------------
check () { [[ -e "$1" ]] || { echo "ERROR: missing $2: $1" >&2; exit 1; }; }
check "${SNPEFF_DIR}/snpEff.jar"                  "snpEff.jar (run setup_resources.sh)"
check "${SNPSIFT_JAR}"                            "SnpSift.jar (run setup_resources.sh)"
check "${SNPEFF_DIR}/snpEff.config"               "snpEff.config"
check "${lofreq_dir}/${name}_snvs.vcf"            "LoFreq SNV file (did step 6 work?)"
check "${SNPEFF_DIR}/data/${SNPEFF_DB}/snpEffectPredictor.bin" "snpEff database ${SNPEFF_DB} (run setup_resources.sh)"
echo "Pre-flight checks OK; $(grep -vc '^#' ${lofreq_dir}/${name}_snvs.vcf || true) SNVs in input"

SNPEFF="java -Xmx8g -jar ${SNPEFF_DIR}/snpEff.jar"
SNPSIFT="java -Xmx4g -jar ${SNPSIFT_JAR}"
CFG="-c ${SNPEFF_DIR}/snpEff.config"

# rCRS region files in the database's chromosome naming (BED, 0-based start)
printf "%s\t16023\t16569\n%s\t0\t576\n" "${SNPEFF_CHROM}" "${SNPEFF_CHROM}" > ${name}_ctrl_region.bed  # 16024-16569 + 1-576
printf "%s\t5720\t5781\n" "${SNPEFF_CHROM}" > ${name}_OriL_region.bed                                  # 5721-5781

FIELDS="CHROM ANN[0].BIOTYPE ANN[0].GENE POS REF ALT QUAL DP DP4[0] DP4[1] DP4[2] DP4[3] AF SB ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P ANN[0].CDS_POS ANN[0].CDS_LEN ANN[0].AA_POS ANN[0].AA_LEN ANN[0].GENEID ANN[0].ERRORS"
HEADER='CHROM\tANN_BIOTYPE\tANN_GENE\tPOS\tREF\tALT\tQUAL\tDP\tDP4_REF_fw\tDP4_REF_rv\tDP4_ALT_fw\tDP4_ALT_rv\tAF\tSB\tANN_EFFECT\tANN_IMPACT\tANN_HGVSc\tANN_HGVSp\tANN_CDSPOS\tANN_CDSLEN\tANN_AAPOS\tANN_AALEN\tANN_GENEID\tANN_ERRORS'

# ---- Filtering ------------------------------------------------------------
# Minimum supporting reads in total and on both strands
$SNPSIFT filter "( (DP*AF >= 15) & (DP4[2] >= 3) & (DP4[3] >= 3) )" \
    ${lofreq_dir}/${name}_snvs.vcf > ${name}_filtered.vcf

# Remove highly strand-biased variants
$SNPSIFT filter "( SB < 1000 )" ${name}_filtered.vcf > ${name}_pos_filtered.vcf

# Rename chrM to the database's mitochondrion name
sed -e "s/^chrM\t/${SNPEFF_CHROM}\t/" -e "s/##contig=<ID=chrM,/##contig=<ID=${SNPEFF_CHROM},/" \
    ${name}_pos_filtered.vcf > ${name}_pos_filtered_${SNPEFF_CHROM}.vcf

# ---- Annotation -------------------------------------------------------------
$SNPEFF ann $CFG -nodownload -no-downstream -no-upstream -noStats ${SNPEFF_DB} \
    ${name}_pos_filtered_${SNPEFF_CHROM}.vcf > ${name}_ann.vcf

# Classic EFF format (gives codon change)
$SNPEFF eff $CFG -nodownload -classic -no-downstream -no-upstream -noStats ${SNPEFF_DB} \
    ${name}_pos_filtered_${SNPEFF_CHROM}.vcf > ${name}_eff.vcf

# ---- Split by genome region -------------------------------------------------
# One annotation per line (replaces snpEff's vcfEffOnePerLine.pl, which is not
# included in recent snpEff downloads): splits the ANN field so each line
# carries a single annotation.
one_per_line () {
    awk 'BEGIN { FS = OFS = "\t" }
    /^#/ { print; next }
    {
        n = split($8, f, ";"); idx = 0
        for (i = 1; i <= n; i++) if (f[i] ~ /^ANN=/) idx = i
        if (idx == 0) { print; next }
        m = split(substr(f[idx], 5), a, ",")
        for (j = 1; j <= m; j++) {
            info = ""
            for (i = 1; i <= n; i++) {
                v = (i == idx) ? "ANN=" a[j] : f[i]
                info = (info == "") ? v : info ";" v
            }
            $8 = info; print
        }
    }'
}

# tRNAs (one effect per line, then keep only tRNA rows)
$SNPSIFT filter "ANN[*].BIOTYPE has 'Mt_tRNA'" ${name}_ann.vcf \
    | one_per_line > ${name}_ann_trna_tmp.vcf

$SNPSIFT filter "ANN[*].BIOTYPE has 'Mt_tRNA'" ${name}_ann_trna_tmp.vcf \
    | $SNPSIFT extractFields - -s "\t" -e "NA" CHROM "ANN[*].BIOTYPE" \
      "ANN[*].GENE" POS REF ALT QUAL DP DP4[0] DP4[1] DP4[2] DP4[3] AF SB "ANN[*].EFFECT" "ANN[*].IMPACT" \
      "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" \
      "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].GENEID" "ANN[*].ERRORS" \
    | sed '1d' | sed "1i${HEADER}" > ${name}_ann_trna.vcf
rm -f ${name}_ann_trna_tmp.vcf

# rRNAs
$SNPSIFT filter "ANN[*].BIOTYPE has 'Mt_rRNA'" ${name}_ann.vcf \
    | $SNPSIFT extractFields - -s "\t" -e "NA" $FIELDS \
    | sed '1d' | sed "1i${HEADER}" > ${name}_ann_rrna.vcf

# Control region
$SNPSIFT intervals ${name}_ctrl_region.bed < ${name}_ann.vcf \
    | $SNPSIFT extractFields - -s "\t" -e "NA" $FIELDS \
    | sed '1d' | sed "1i${HEADER}" > ${name}_ann_ctrl.vcf

# OriL
$SNPSIFT intervals ${name}_OriL_region.bed < ${name}_ann.vcf \
    | $SNPSIFT extractFields - -s "\t" -e "NA" $FIELDS \
    | sed '1d' | sed "1i${HEADER}" > ${name}_ann_oril.vcf

# Protein coding
$SNPSIFT filter "ANN[*].BIOTYPE has 'protein_coding'" ${name}_ann.vcf \
    | $SNPSIFT extractFields - -s "\t" -e "NA" $FIELDS \
    | sed '1d' | sed "1i${HEADER}" > ${name}_ann_protein_tmp.vcf

# Codon position (1st/2nd/3rd) from ANN_CDSPOS (column 19)
sed '1d' ${name}_ann_protein_tmp.vcf | cut -f 19 | while read -r i; do
    if [[ -z "$i" || ! "$i" =~ ^[0-9]+$ ]]; then echo "NA"
    elif (( i % 3 == 1 )); then echo "1"
    elif (( i % 3 == 2 )); then echo "2"
    else echo "3"
    fi
done | sed '1iCODON_POS' > ${name}_codon_position.tmp

# Codons and amino-acid change from the classic EFF annotation
$SNPSIFT filter "EFF[*].BIOTYPE has 'protein_coding'" ${name}_eff.vcf \
    | $SNPSIFT extractFields - -s "\t" -e "NA" "EFF[0].CODON" "EFF[0].AA" \
    | sed '/^#/d;s/\//\t/g' \
    | sed '1iEFF_CODONREF\tEFF_CODONALT\tEFF_HGVSp' > ${name}_protein_eff.tmp

paste ${name}_ann_protein_tmp.vcf ${name}_codon_position.tmp ${name}_protein_eff.tmp > ${name}_ann_protein.vcf
rm -f ${name}_ann_protein_tmp.vcf ${name}_codon_position.tmp ${name}_protein_eff.tmp
rm -f ${name}_ctrl_region.bed ${name}_OriL_region.bed

# ---- Combine all regions into one genome-wide table ------------------------
sed -e '2,${/^CHROM/d' -e '}' \
    ${name}_ann_ctrl.vcf ${name}_ann_oril.vcf ${name}_ann_protein.vcf \
    ${name}_ann_rrna.vcf ${name}_ann_trna.vcf \
    | sed '1s/$/\tCODON_POS\tEFF_CODONREF\tEFF_CODONALT\tEFF_HGVSp/' > ${name}_ann_MT.vcf

echo "Annotation done for ${name}"
