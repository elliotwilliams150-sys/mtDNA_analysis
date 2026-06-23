#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(ape)
  library(phylobase)
  library(phylosignal)
  library(RhpcBLASctl)
  library(parallel)
  library(dplyr)
  library(tibble)
  library(data.table)
  library(glue)
})

set.seed(123)
repo_dir = '/lab-share/Hem-Sankaran-e2/Public/projects/tgao/mitodrift_analysis'
source(glue("{repo_dir}/R/lineage_coupling.R"))

#
# Argument parsing
#
parser <- ArgumentParser(description = "Run phylogenetic signal analysis on trait matrix")
parser$add_argument("--trait_mat", required=TRUE, help="Path to trait matrix file (CSV) - traits in rows, cells in columns")
parser$add_argument("--tree", required=TRUE, help="Path to phylogenetic tree file (Newick format)")
parser$add_argument("--outfile", required=TRUE, help="Output file path for results")
parser$add_argument("--reps", type="integer", default=999, help="Number of permutations for phyloSignal [default: %(default)s]")
parser$add_argument("--ncores", type="integer", default=1, help="Number of cores to use [default: %(default)s]")

args <- parser$parse_args()

# Prepare output folders
base_dir <- dirname(args$outfile)
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive=TRUE, showWarnings=FALSE)
}

# Load trait matrix
message("Loading trait matrix...")
trait_mat <- fread(args$trait_mat) %>% tibble::column_to_rownames('V1') %>% as.matrix
message(sprintf("Loaded trait matrix: %d traits x %d cells", nrow(trait_mat), ncol(trait_mat)))

# Load phylogenetic tree
message("Loading phylogenetic tree...")
tree <- read.tree(args$tree)
message(sprintf("Loaded tree with %d tips", length(tree$tip.label)))

# Take intersection of cells in phylo and trait mat
message("Taking intersection of cells...")
cells_common <- intersect(colnames(trait_mat), tree$tip.label)
message(sprintf("Found %d common cells between trait matrix and tree", length(cells_common)))

# Filter trait matrix to common cells
trait_mat <- trait_mat[, cells_common, drop = FALSE]
message(sprintf("Filtered trait matrix to common cells: %d traits x %d cells", nrow(trait_mat), ncol(trait_mat)))

# Filter tree to common cells
tree <- tree %>% ape::keep.tip(cells_common)
tree$node.label <- NULL
message(sprintf("Filtered tree to common cells: %d tips", length(tree$tip.label)))

# Filter out features with zero variance among cells
message("Filtering out features with zero variance...")
feature_variances <- apply(trait_mat, 1, var, na.rm = TRUE)
zero_var_features <- which(feature_variances == 0)
if (length(zero_var_features) > 0) {
  message(sprintf("Removing %d features with zero variance", length(zero_var_features)))
  trait_mat <- trait_mat[-zero_var_features, , drop = FALSE]
  message(sprintf("After variance filtering: %d traits x %d cells", nrow(trait_mat), ncol(trait_mat)))
} else {
  message("No features with zero variance found")
}

# Run phylogenetic signal analysis
message("Running phylogenetic signal analysis...")
results <- run_signal(
  trait_mat = trait_mat,
  tree = tree,
  reps = args$reps,
  ncores = args$ncores
)

# Write results
fwrite(results, args$outfile)
message(sprintf("Wrote results to: %s", args$outfile))

# Print summary
message(sprintf("Analysis complete. Found %d traits with significant phylogenetic signal (qval < 0.05)", 
                sum(results$qval < 0.05, na.rm = TRUE)))
