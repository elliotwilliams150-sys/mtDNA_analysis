#!/bin/bash
#random sampler, mut_dat generator 
#N_CELLS=200

#SUBSAMPLE_R=$(mktemp /tmp/subsample_XXXX.R)
#cat > "$SUBSAMPLE_R" << 'EOF'
#args <- commandArgs(trailingOnly = TRUE)
#dat <- read.csv(args[1], header = FALSE, col.names = c("cell", "variant", "a", "d"))
#set.seed(as.integer(args[4]))
#cells_keep <- sample(unique(dat$cell), as.integer(args[3]))
#sub <- dat[dat$cell %in% cells_keep, ]
#write.csv(sub, args[2], row.names = FALSE)
#cat("Wrote", nrow(sub), "rows for", length(cells_keep), "cells\n")
#EOF

#YIPEEEEE

echo "=== PMMono ==="

for SEED in $(seq 1 20); do
    OUTDIR="/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets/PMMono/subsample_seed${SEED}"
    mkdir -p "$OUTDIR"
    SUB_CSV="${OUTDIR}/mut_dat_sub.csv"
    echo "Seed $SEED"
    # optional if subsample script ran Rscript "$SUBSAMPLE_R" "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets/PMMono" "$SUB_CSV" "$N_CELLS" "$SEED"
    Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
        --mut_dat "$SUB_CSV" --outdir "$OUTDIR" \
        --tree_mcmc_iter 200000 --tree_mcmc_chains 50 --tree_mcmc_burnin 10000 \
        --ncores 50 --ncores_em 20 --ncores_nj 20 --ncores_qs 20 --ncores_annot 1 --resume TRUE
done

echo "=== SurgMono ==="
for SEED in $(seq 1 20); do
    OUTDIR="/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets/SurgMono/subsample_seed${SEED}"
    mkdir -p "$OUTDIR"
    SUB_CSV="${OUTDIR}/mut_dat_sub.csv"
    echo "Seed $SEED"
    # optional if subsample script ran Rscript "$SUBSAMPLE_R" "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets/SurgMono" "$SUB_CSV" "$N_CELLS" "$SEED"
    Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
        --mut_dat "$SUB_CSV" --outdir "$OUTDIR" \
        --tree_mcmc_iter 100000 --tree_mcmc_chains 50 --tree_mcmc_burnin 10000 \
        --ncores 50 --ncores_em 20 --ncores_nj 20 --ncores_qs 20 --ncores_annot 1 --resume FALSE
done

rm "$SUBSAMPLE_R"






