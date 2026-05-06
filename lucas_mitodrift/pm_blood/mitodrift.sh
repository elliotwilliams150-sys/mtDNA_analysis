module load GLPK
module load R

Rscript ../../../Programs/run_mitodrift_em.R \
  --mut_dat PBMC_PM_long_variant_matrix.csv \
  --outdir PBMC_PM \
  --tree_mcmc_iter 2000 \
  --tree_mcmc_chains 4 \
  --tree_mcmc_burnin 200 \
  --ncores 8 
