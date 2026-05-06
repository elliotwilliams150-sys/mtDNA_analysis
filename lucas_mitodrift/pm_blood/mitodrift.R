library(mitodrift)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtree)

# 1. Your mutation data (same one you used as input)
mut_dat <- read.csv("PBMC_PM_long_variant_matrix.csv")

# 2. Your model output
md <- readRDS("/nobackup/proj/rockhpc_spnmmd/LucasCortes/Mike_keogh_mtscATAC/mitodrift/pm_blood/PBMC_PM/mitodrift_object.rds")

# 3. Trim the tree
phy_trim <- trim_tree(md$tree_annot, conf = 0.5)

# 4. Plot
pdf("PBMC_PM_Heatmap1.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  phy_trim,
  mut_dat,
  dot_size = 1,
  branch_width = 0.3,
  branch_length = FALSE,
  node_conf = TRUE,
  het_max = 1,
  title = "PBMC_postmortem_tree"
)
dev.off()
pdf("precision_recall_PBMC_PM.pdf")
pr_df <- compute_variant_pr_curve(md$tree_annot, mut_dat)
plot_prec_recall_vs_conf(
  pr_df,
  sample_name = "Variant-based precision recall PBMC_PM",
  cutoff = 0.3
)
dev.off()
# 5. Trim tree based on the confidence threshold of previous plot

tree_trim <- trim_tree(md$tree_annot, conf = 0.3)

pdf("Trimmed_tree_0.3_PM_PBMC.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  node_conf = TRUE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Trimmed tree (conf >= 0.3) PBMC_PM"
)
dev.off()
# 6. Save the variants in the trimmed tree
tree_variants = unique(mut_dat$variant)
write.csv(tree_variants, "tree_variants_PBMC_PM.csv", row.names = FALSE)
#head(mut_dat, n=20)
# 7. Assign clones

clone_df <- assign_clones_polytomy(tree_trim)


# 8. Visualise the clones
clade_order <- unique(clone_df$clade)
clone_pal <- make_clade_pal(length(clade_order), labels = clade_order,
                            pal = "Dark2", cycle_len = 8, cycle_shift = 0)
pdf("Clone_tree_0.3_PBMC_PM.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  cell_annot = clone_df,
  annot_pal = clone_pal,
  node_conf = TRUE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Clones on trimmed tree PBMC_PM"
)
dev.off()
head(clone_df)
head(clone_pal)
pdf("Circ_chart_0.3_PBMC_PM.pdf", width = 10, height = 8)
plot_phylo_circ(
  tree_trim,
  cell_annot = clone_df,
  annot_pal = clone_pal,
  annot_legend = FALSE,
  title = "Circular layout"
)
dev.off()
