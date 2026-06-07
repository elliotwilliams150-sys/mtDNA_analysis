#!/bin/bash
N_CELLS=200

SUBSAMPLE_R=$(mktemp /tmp/subsample_XXXX.R)
cat > "$SUBSAMPLE_R" << 'EOF'
args <- commandArgs(trailingOnly = TRUE)
dat <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE)
colnames(dat) <- c("cell", "variant", "a", "d")  # Standardize names
dat$a <- as.numeric(dat$a)
dat$d <- as.numeric(dat$d)
dat <- dat[complete.cases(dat[, c("a", "d")]) & dat$d > 0 & dat$a >= 0, ]
set.seed(as.integer(args[4]))
n_keep <- min(as.integer(args[3]), length(unique(dat$cell)))
cells_keep <- sample(unique(dat$cell), n_keep)
sub <- dat[dat$cell %in% cells_keep, ]
write.csv(sub, args[2], row.names = FALSE, quote = FALSE)
cat("Wrote", nrow(sub), "rows for", length(cells_keep), "cells from", nrow(dat), "total\n")
EOF

echo "=== SurgINFLAM ==="
for SEED in $(seq 1 20); do
    OUTDIR="/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes/INFLAM/subsample_seed${SEED}"
    mkdir -p "$OUTDIR"
    SUB_CSV="${OUTDIR}/mut_dat_sub.csv"
    echo "Seed $SEED"
    Rscript "$SUBSAMPLE_R" "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes/INFLAM/INFLAM.csv" "$SUB_CSV" "$N_CELLS" "$SEED"
    Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
        --mut_dat "$SUB_CSV" --outdir "$OUTDIR" \
        --tree_mcmc_iter 35000 --tree_mcmc_chains 20 --tree_mcmc_burnin 5000 \
        --ncores 20 --ncores_em 20 --ncores_nj 20 --ncores_qs 20 --ncores_annot 1 --resume FALSE
done

echo "=== SurgHOMEO ==="
for SEED in $(seq 1 20); do
    OUTDIR="/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes/HOMEO/subsample_seed${SEED}"
    mkdir -p "$OUTDIR"
    SUB_CSV="${OUTDIR}/mut_dat_sub.csv"
    echo "Seed $SEED"
    Rscript "$SUBSAMPLE_R" "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes/HOMEO/HOMEO.csv" "$SUB_CSV" "$N_CELLS" "$SEED"
    Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
        --mut_dat "$SUB_CSV" --outdir "$OUTDIR" \
        --tree_mcmc_iter 35000 --tree_mcmc_chains 20 --tree_mcmc_burnin 5000 \
        --ncores 20 --ncores_em 20 --ncores_nj 20 --ncores_qs 20 --ncores_annot 1 --resume FALSE
done

rm "$SUBSAMPLE_R"