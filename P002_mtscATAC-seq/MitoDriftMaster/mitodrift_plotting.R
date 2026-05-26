#whoopee
library(mitodrift)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtree)

# 1. Your mutation data (same one you used as input)
mut_dat <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/combined_mut_dat.csv")

# 2. Your model output
md <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/mitodrift_object.rds")

# 3. Optional - tree diangostics 
diag <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/tree_mcmc_diag.rds")
str(diag) 
print(diag$asdsf)


# 3. Trim the tree
tree_trim <- trim_tree(md$tree_annot, conf = 0.125)

# 4. Plot
pdf("MonoSubset1.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  md$tree_annot,
  mut_dat,
  dot_size = 1,
  branch_width = 0.3,
  branch_length = FALSE,
  node_conf = FALSE,
  node_lab = FALSE,
  conf_label = TRUE,  
  het_max = 1,
  ladderize =  TRUE, 

  title = "MonoSubset1_tree"
)
dev.off()

pr_df <- compute_variant_pr_curve(md$tree_annot, mut_dat)
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/PR_curve.pdf", width = 6, height = 4)
plot_prec_recall_vs_conf(
  pr_df,
  sample_name = "Variant-based precision recall Ctx_Surgical",
  cutoff = 0.18
)
dev.off()
# 5. Trim tree based on the confidence threshold of previous plot

tree_trim <- trim_tree(md$tree_annot, conf = 0.18)

pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/Trimmed_tree_Ctx_Surgical.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  node_conf = TRUE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Trimmed tree MonoSubset1"
)
dev.off()
# 6. Save the variants in the trimmed tree
tree_variants = unique(mut_dat$variant)
write.csv(tree_variants, "tree_variants_Ctx_Surgical.csv", row.names = FALSE)
#head(mut_dat, n=20)
# 7. Assign clones

clone_df <- assign_clones_polytomy(tree_trim)

# 8. Visualise the clones
clade_order <- unique(clone_df$clade)
clone_pal <- make_clade_pal(length(clade_order), labels = clade_order,
                            pal = "Dark2", cycle_len = 8, cycle_shift = 0)
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/Clone_tree_Ctx_Surgical.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  cell_annot = list(clone_df, df_with_annot, cell_type_df), 
  annot_pal = list(clone_pal, c(P = "gray50", H = "orange"), c(P = "white", H = "white", NK = "white", Mono = "white", B_cell = "white", T_cell = "white",
  Oligos = "white",
  HOMEO = "green",
  INFLAM = "red",
  `T cell` = "white"
)   ),
  layered = FALSE,
  node_conf = FALSE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Clones on trimmed tree",
)
dev.off()
head(clone_df)
head(clone_pal)
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets/PMMono/subsample_seed1/Circ_chart_0.1_Ctx_Surgical.pdf", width = 10, height = 8)
plot_phylo_circ(
  tree_trim,
  cell_annot = clone_df,
  annot_pal = clone_pal,
  annot_legend = FALSE,
  title = "Circular layout"
)
dev.off()



df_with_annot <- clone_df %>%
  mutate(
    annot = substr(gsub("^(.)(_.*)", "\\1", cell), 1, 1),
    annot = if_else(annot %in% c("P", "H"), annot, "Unknown")  # safety if you ever see non‑P/H
  )

df_with_annot

cell_type_df <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/tissue_mitodrift_df_combined.csv", stringsAsFactors = FALSE)

#drop singletons 
drop_root_singletons <- function(tree) {
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
  root_children <- tree$edge[tree$edge[, 1] == root, 2]
  Ntip <- length(tree$tip.label)
  singleton_tips <- root_children[root_children <= Ntip]
  
  if (length(singleton_tips) == 0) return(tree)
  
  ape::drop.tip(tree, tree$tip.label[singleton_tips])
}

tree_trim <- drop_root_singletons(tree_trim)

#Clade enrichment analysis 
library(ape)
library(phangorn)

