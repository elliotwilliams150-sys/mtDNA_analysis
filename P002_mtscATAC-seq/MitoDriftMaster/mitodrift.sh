  #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDriftMaster/Combined_Surgical_Lineage_Matrix.csv \
  --outdir /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526 \
  --tree_mcmc_iter 20000 \
  --tree_mcmc_chains 10 \
  --tree_mcmc_burnin 2000 \
  --ncores 10 \
  --ncores_em 10 \
  --ncores_nj 10 \
  --ncores_qs 10 \
  --ncores_annot 1 \
  --resume FALSE






