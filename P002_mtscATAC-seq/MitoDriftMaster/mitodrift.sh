 #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526/combined_mut_dat.csv \
  --outdir  /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526 \
  --tree_mcmc_iter 250000 \
  --tree_mcmc_chains 10 \
  --tree_mcmc_burnin 20000 \
  --conv_thres 0.1 \
  --ncores 1 \
  --ncores_em 1 \
  --ncores_nj 1 \
  --ncores_qs 1 \
  --ncores_annot 1 \
  --resume TRUE
  