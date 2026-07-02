 #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/PM_combined_matrix.csv \
  --outdir  /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/ \
  --tree_mcmc_iter 500000 \
  --tree_mcmc_chains 10 \
  --tree_mcmc_burnin 50000 \
  --conv_thres 0.05 \
  --ncores 10 \
  --ncores_em 10 \
  --ncores_nj 10 \
  --ncores_qs 10 \
  --ncores_annot 10 \
  --resume FALSE
  