 #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets/PM/PM_Monos_subset.csv \
  --outdir /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/PM \
  --tree_mcmc_iter 35000 \
  --tree_mcmc_chains 30 \
  --tree_mcmc_burnin 3500 \
  --ncores 30 \
  --ncores_em 30 \
  --ncores_nj 30 \
  --ncores_qs 30 \
  --ncores_annot 1 \
  --resume FALSE