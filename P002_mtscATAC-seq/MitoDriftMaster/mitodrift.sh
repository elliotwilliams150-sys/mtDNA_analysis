 #!/usr/bin/env bash
  #module load GLPK  # Loads libglpk for optimization
  #module load R

cd /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526

Rscript /home/elliot/R/x86_64-pc-linux-gnu-library/4.4/mitodrift/bin/run_mitodrift_em.R \
  --mut_dat /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Surg_McgVSAst_mtSDI_sub/AST/AST.csv \
  --outdir /mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Surg_McgVSAst_mtSDI_sub/AST \
  --tree_mcmc_iter 10000 \
  --tree_mcmc_chains 20 \
  --tree_mcmc_burnin 1000 \
  --ncores 20 \
  --ncores_em 30 \
  --ncores_nj 30 \
  --ncores_qs 30 \
  --ncores_annot 1 \
  --resume FALSE