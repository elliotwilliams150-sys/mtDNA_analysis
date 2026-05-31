 #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Cluser7Cat

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Cluser7Cat/combined_mut_dat.csv \
  --outdir  /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Cluser7Cat \
  --tree_mcmc_iter 200000  \
  --tree_mcmc_chains 10 \
  --tree_mcmc_burnin 25000 \
  --conv_thres 0.1296 \
  --ncores 2 \
  --ncores_em 10 \
  --ncores_nj 10 \
  --ncores_qs 1 \
  --ncores_annot 1 \
  --resume TRUE
  