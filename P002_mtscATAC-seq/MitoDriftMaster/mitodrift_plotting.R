#whoopee
library(mitodrift)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtree)

# 1. Your mutation data (same one you used as input)
mut_dat <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/Combined_Surgical_Lineage_Matrix.csv")

# 2. Your model output
md <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/mitodrift_object.rds")

# 3. Optional - tree diangostics 
diag <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/tree_mcmc_diag.rds")
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
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/PR_curve.pdf", width = 6, height = 4)
plot_prec_recall_vs_conf(
  pr_df,
  sample_name = "Variant-based precision recall Ctx_Surgical",
  cutoff = 0.60
)
dev.off()
# 5. Trim tree based on the confidence threshold of previous plot

tree_trim <- trim_tree(md$tree_annot, conf = 0.85)

pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/Trimmed_tree_Ctx_Surgical.pdf", width = 10, height = 8)
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
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/Clone_tree_Ctx_Surgical.pdf", width = 10, height = 8)
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


#drop singletons 
tree_trim <- drop_root_singletons(tree_trim)

drop_root_singletons <- function(tree) {
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
  root_children <- tree$edge[tree$edge[, 1] == root, 2]
  Ntip <- length(tree$tip.label)
  singleton_tips <- root_children[root_children <= Ntip]
  
  if (length(singleton_tips) == 0) return(tree)
  
  ape::drop.tip(tree, tree$tip.label[singleton_tips])
}

cell_type_df <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun040526/tissue_mitodrift_df_combined.csv", stringsAsFactors = FALSE)

#Clade enrichment analysis 
trait <- cell_type_df$annot
names(trait) <- cell_type_df$cell
trait <- trait[tree_trim$tip.label]
anyNA(trait)
table(trait)

trait_bin <- as.integer(trait %in% c("Mono"))
names(trait_bin) <- names(trait)

global_mean <- mean(trait_bin, na.rm = TRUE)
global_sd <- sd(trait_bin, na.rm = TRUE)

nodes <- (length(tree_trim$tip.label) + 1):(length(tree_trim$tip.label) + tree_trim$Nnode)

res_mono <- do.call(rbind, lapply(nodes, function(node) {
  desc <- phangorn::Descendants(tree_trim, node, type = "tips")[[1]]
  tip_names <- tree_trim$tip.label[desc]
  n <- length(tip_names)
  if (n < 5) return(NULL)

  clade_mean <- mean(trait_bin[tip_names], na.rm = TRUE)
  deviation <- clade_mean - global_mean
  se <- global_sd / sqrt(n)

  if (is.na(se) || se == 0) {
    z <- NA_real_
    p <- NA_real_
  } else {
    z <- deviation / se
    p <- 2 * pnorm(-abs(z))
  }

  data.frame(
    node = node,
    size = n,
    clade_mean = clade_mean,
    deviation = deviation,
    z = z,
    p = p
  )
}))

res_mono$padj <- p.adjust(res_mono$p, method = "BH")

head(res_mono[order(-res_mono$z), ])
head(res_mono[order(res_mono$z), ])
summary(res_mono$padj)


#we have already built a binary trait vector, global mean and SD,  per-node clade means, deviations, z-scores, p-values, and BH-adjusted p-values
cmp <- merge(
  merge(res_mono[, c("node", "z")], res_inflam[, c("node", "z")],
        by = "node", suffixes = c("_mono", "_inflam")),
  res_homeo[, c("node", "z")],
  by = "node"
)
names(cmp)[4] <- "z_homeo"

cor(cmp$z_mono, cmp$z_inflam, use = "complete.obs")
cor(cmp$z_mono, cmp$z_homeo, use = "complete.obs")

cmp$diff_inflam <- abs(cmp$z_mono - cmp$z_inflam)
cmp$diff_homeo  <- abs(cmp$z_mono - cmp$z_homeo)

mean(cmp$diff_inflam, na.rm = TRUE)
mean(cmp$diff_homeo, na.rm = TRUE)

median(cmp$diff_inflam, na.rm = TRUE)
median(cmp$diff_homeo, na.rm = TRUE)

#inter-phylogneetic differences 
library(picante)
library(ape)
library(picante)
library(dplyr)

# -------------------------------
# Inputs:
# tree_trim    : phylo object
# cell_type_df : data.frame with columns cell and annot
# -------------------------------

# Map annotations onto tree tips
trait <- cell_type_df$annot
names(trait) <- cell_type_df$cell
trait <- trait[tree_trim$tip.label]

anyNA(trait)
table(trait)

# Set all branch lengths to 1, as in the paper
tree_unit <- tree_trim
tree_unit$edge.length <- rep(1, nrow(tree_unit$edge))

# Pairwise tip-to-tip distances
phy_dist <- cophenetic(tree_unit)

# Groups to compare
groups <- c("Mono", "INFLAM", "HOMEO")
groups <- groups[groups %in% unique(trait)]

# Group-by-tip incidence matrix
comm <- sapply(groups, function(g) as.integer(trait == g))
comm <- t(comm)
colnames(comm) <- names(trait)

# Match tip order to phylogenetic distance matrix
comm <- comm[, colnames(phy_dist), drop = FALSE]

# Inter-community phylogenetic distances
cd <- picante::comdist(comm, phy_dist)

cd

# Direct comparison for Mono
if ("Mono" %in% rownames(cd)) {
  cd["Mono", setdiff(colnames(cd), "Mono")]
}

# Optional direct printout if all three groups are present
if (all(c("Mono", "INFLAM", "HOMEO") %in% rownames(cd))) {
  cat("Mono vs INFLAM:", cd["Mono", "INFLAM"], "\n")
  cat("Mono vs HOMEO :", cd["Mono", "HOMEO"], "\n")
}