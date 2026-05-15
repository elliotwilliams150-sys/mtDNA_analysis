 #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun150526/Combined_Surgical_Lineage_Matrix.csv \
  --outdir  /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun150526 \
  --tree_mcmc_iter 100000 \
  --tree_mcmc_chains 10 \
  --tree_mcmc_burnin 10000 \
  --conv_thresh 0.1 \
  --ncores 10 \
  --ncores_em 30 \
  --ncores_nj 30 \
  --ncores_qs 30 \
  --ncores_annot 1 \
  --resume FALSE