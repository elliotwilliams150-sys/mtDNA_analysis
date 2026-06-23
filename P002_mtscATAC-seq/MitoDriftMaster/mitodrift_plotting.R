#whoopee
library(mitodrift)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtree)
source("/mnt/claw-raid/elliot/P002_mtscATAC-seq/scripts/path_state_fate_analysis/R/utils.R")
source("/mnt/claw-raid/elliot/P002_mtscATAC-seq/scripts/path_state_fate_analysis/R/stats.R")
source("/mnt/claw-raid/elliot/P002_mtscATAC-seq/scripts/path_state_fate_analysis/R/plotting_themes.R")
source("/mnt/claw-raid/elliot/P002_mtscATAC-seq/scripts/path_state_fate_analysis/R/lineage_coupling.R")

# 1. Your mutation data (same one you used as input)
mut_dat <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/PM_combined_matrix.csv")

# 2. Your model output
md <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/mitodrift_object.rds")

# 3. Optional - tree diangostics 
diag <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRunAD_copy/tree_mcmc_diag.rds")
str(diag) 
print(diag$asdsf)


# 3. Trim the tree
tree_trim <- trim_tree(md$tree_annot, conf = 0.05)

# 4. Plot
pdf("PMcat.pdf", width = 10, height = 8)
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

  title = "PM_tree"
)
dev.off()

pr_df <- compute_variant_pr_curve(md$tree_annot, mut_dat)
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/PR_curve.pdf", width = 6, height = 4)
plot_prec_recall_vs_conf(
  pr_df,
  sample_name = "Variant-based precision recall Ctx_Surgical",
  cutoff = 0.03
)
dev.off()
# 5. Trim tree based on the confidence threshold of previous plot

tree_trim <- trim_tree(md$tree_annot, conf = 0.15)


pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/Trimmed_tree.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  node_conf = TRUE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Trimmed tree Blood+Brain_PM"
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
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRunAD_copy/Clone_tree_BainBlood_AD.pdf", width = 10, height = 8)

plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  cell_annot = list(clone_df, cell_type_df),
  annot_pal = list(
    clone_pal,
    c(
      Monocytes = "#ffffff",
      B_cell = "white",
      T_cell = "white",
      Oligo = "white",
      Microglia = "#ffffff",
      `Marco+` = "red",
      `Mcg/PVM` = "#ffffff"
    )
  ),
  layered = FALSE,
  node_conf = FALSE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Clones on trimmed tree"
)

dev.off()
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


pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRunAD_copy/Clone_tree_BainBlood_AD.pdf", width = 10, height = 8)

plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  cell_annot = list(clone_df, cell_type_df),
  annot_pal = list(
    clone_pal,
    annot_colors
  ),
  layered = FALSE,
  node_conf = FALSE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Clones on trimmed tree"
)

dev.off()

pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/ForHPC/PM/Clone_tree.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  cell_annot = list(clone_df, df_with_annot), 
  annot_pal = list(clone_pal, c(P = "gray50", C = "orange")),
  layered = FALSE,
  node_conf = FALSE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Clones on trimmed tree",
)
dev.off()

df_with_annot <- clone_df %>%
  mutate(
    annot = substr(gsub("^(.)(_.*)", "\\1", cell), 1, 1),
    annot = if_else(annot %in% c("P", "C"), annot, "Unknown")  # safety if you ever see non‑P/H
  )

df_with_annot


annot_colors <- c(
  Oligo = "white",
  T_cell = "white",
  Microglia = "white",
  `Mcg/PVM` = "white",
  `Marco+` = "red",
  B_cell = "white",
  CD8 = "white",
  CD4 = "white",
  Monocytes = "blue"
)


#create cell type df
# load the Leiden cluster files (must have cell + cluster column)
brain_clusters <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRunAD_copy/brain_cell_leiden_clusters.csv",
                           stringsAsFactors = FALSE)
blood_clusters <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRunAD_copy/blood_cell_leiden_clusters.csv",
                           stringsAsFactors = FALSE)

# add the prefix so they match clone_df
brain_clusters$cell <- paste0("H_", brain_clusters$cell)
blood_clusters$cell <- paste0("P_", blood_clusters$cell)

# add tissue labels
brain_clusters$tissue <- "brain"
blood_clusters$tissue <- "blood"

# combine
leiden_cluster_df <- rbind(brain_clusters, blood_clusters)

# make annot include tissue prefix
leiden_cluster_df$annot <- ifelse(
  leiden_cluster_df$tissue == "brain",
  paste0("H_", leiden_cluster_df$leiden),
  paste0("P_", leiden_cluster_df$leiden)
)

# keep cell as it is, with prefix already in it
# if needed, make sure cell has the right prefixes too
# brain should be H_..., blood should be P_...

# inspect
colnames(leiden_cluster_df)
head(leiden_cluster_df)

# build final dataframe
cell_type_df <- merge(
  clone_df[, c("cell", "clade", "clade_node", "size", "frac")],
  leiden_cluster_df[, c("cell", "annot", "tissue")],
  by = "cell",
  all.x = TRUE
)

cell_type_df <- cell_type_df[, c("cell", "clade", "clade_node", "annot", "size", "frac", "tissue")]

cell_type_df$annot <- as.character(cell_type_df$annot)

cell_type_df$annot[cell_type_df$annot == "H_0"] <- "Oligo"
cell_type_df$annot[cell_type_df$annot %in% c("H_1", "H_3")] <- "T cell"
cell_type_df$annot[cell_type_df$annot == "H_4"] <- "Microglia"
cell_type_df$annot[cell_type_df$annot == "H_2"] <- "P Mcg/PVM?"
cell_type_df$annot[cell_type_df$annot == "H_5"] <- "Marco+"
cell_type_df$annot[cell_type_df$annot == "P_6"] <- "B cell"
cell_type_df$annot[cell_type_df$annot == "P_4"] <- "CD8"
cell_type_df$annot[cell_type_df$annot == "P_7"] <- "CD4"
cell_type_df$annot[cell_type_df$annot %in% c("P_1", "P_2", "P_3")] <- "Monocytes"


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

#Create tissue type annotation dataframe
df_with_annot <- mut_dat %>%
  mutate(
    annot = substr(gsub("^(.)(_.*)", "\\1", cell), 1, 1),
    annot = if_else(annot %in% c("P", "C"), annot, "Unknown")  # safety if you ever see non‑P/H
  ) 

# build cell_type_df in the format you want
cell_type_df <- merge(
  clone_df[, c("cell", "clade", "clade_node", "size", "frac")],
  leiden_cluster_df[, c("cell", "leiden", "tissue")],
  by = "cell",
  all.x = TRUE
)

cell_type_df <- cell_type_df[, c("cell", "clade", "clade_node", "leiden", "size", "frac", "tissue")]



#Create cell type annotation dataframe
cell_type_df <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/tissue_mitodrift_df_combined.csv", stringsAsFactors = FALSE)





                                                    ### clade enrichment analysis ### 
shared <- intersect(tree_trim$tip.label, cell_type_df$cell)
tree_trim2 <- ape::keep.tip(tree_trim, shared)
cell_type_df2 <- cell_type_df[cell_type_df$cell %in% shared, ]
trait <- setNames(cell_type_df2$annot, cell_type_df2$cell)
trait <- trait[tree_trim2$tip.label]
mono <- as.numeric(trait == "INFLAM")
names(mono) <- names(trait)
# binary monocyte trait: 1 = Mono, 0 = not Mono
mono <- as.numeric(trait == "INFLAM")
names(mono) <- tree_trim2$tip.label
# global mean and SD across all tips
global_mean <- mean(mono, na.rm = TRUE)
global_sd <- sd(mono, na.rm = TRUE)
# get internal nodes
Ntip <- length(tree_trim2$tip.label)
nodes <- (Ntip + 1):(Ntip + tree_trim2$Nnode)
# compute clade mean and deviation for each internal node
mono_clade_stats <- data.frame(
  node = nodes,
  size = NA_integer_,
  clade_mean = NA_real_,
  deviation = NA_real_
)
for (i in seq_along(nodes)) {
  node <- nodes[i]
  desc <- phangorn::Descendants(tree_trim2, node, type = "tips")[[1]]
  tip_names <- tree_trim2$tip.label[desc]
  vals <- mono[tip_names]
  
  mono_clade_stats$size[i] <- length(vals)
  mono_clade_stats$clade_mean[i] <- mean(vals, na.rm = TRUE)
  mono_clade_stats$deviation[i] <- mono_clade_stats$clade_mean[i] - global_mean
}
head(mono_clade_stats)

#Deviation from gloabl mean / size-adjusted SE = z score, use clone_df to get clade size
mono_clade_stats <- clone_df %>%
  group_by(clade) %>%
  summarise(
    size = n(),
    clade_mean = mean(mono[cell], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    deviation = clade_mean - global_mean,
    se = global_sd / sqrt(size),
    z = deviation / se
  )
head(mono_clade_stats)

#calculate p value from z score 
min_clade_size <- 10 #select major clades, do this per cell type
mono_clade_stats$p <- NA_real_
keep <- !is.na(mono_clade_stats$z) & mono_clade_stats$size >= min_clade_size
# two-sided p-values from standard normal
mono_clade_stats$p[keep] <- 2 * pnorm(-abs(mono_clade_stats$z[keep]))
# Benjamini-Hochberg correction
mono_clade_stats$padj <- NA_real_
mono_clade_stats$padj[keep] <- p.adjust(mono_clade_stats$p[keep], method = "BH")
head(mono_clade_stats)


#'clade enrichment'
library(ape)
library(dplyr)

# -----------------------------
# 1) Prepare trimmed tree and matching cell annotation
# -----------------------------
shared <- intersect(tree_trim$tip.label, cell_type_df$cell)
tree_trim2 <- ape::keep.tip(tree_trim, shared)
cell_type_df2 <- cell_type_df[cell_type_df$cell %in% shared, ]

# -----------------------------
# 2) Use existing mono_clade_stats, keep only positive enrichment
#    Expected columns: clade, z, size, padj
# -----------------------------
mono_clade_stats2 <- mono_clade_stats %>%
  distinct(clade, .keep_all = TRUE) %>%
  mutate(clade = as.character(clade)) %>%
  filter(size >= 10, z > 0)
# Optional stricter version:
# mono_clade_stats2 <- mono_clade_stats %>%
#   distinct(clade, .keep_all = TRUE) %>%
#   mutate(clade = as.character(clade)) %>%
#   filter(size >= 10, z > 0, !is.na(padj), padj < 0.05)

# -----------------------------
# 3) Map each clade to its MRCA node in the trimmed tree
# -----------------------------
clade_tips <- clone_df %>%
  group_by(clade) %>%
  summarise(
    tips = list(intersect(unique(cell), tree_trim2$tip.label)),
    .groups = "drop"
  )

clade_node_map2 <- lapply(seq_len(nrow(clade_tips)), function(i) {
  tips_i <- clade_tips$tips[[i]]
  if (length(tips_i) < 2) return(NULL)

  node_i <- ape::getMRCA(tree_trim2, tips_i)
  if (is.null(node_i) || is.na(node_i)) return(NULL)

  data.frame(
    clade = as.character(clade_tips$clade[i]),
    clade_node = as.integer(node_i),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

tmp <- merge(
  clade_node_map2,
  mono_clade_stats2[, c("clade", "z", "size", "padj")],
  by = "clade",
  all.x = TRUE
)

tmp <- tmp[!is.na(tmp$z), ]
tmp <- tmp[tmp$size >= 10, ]

cat("Positive enriched clades mapped:", nrow(tmp), "\n")

# -----------------------------
# 4) Recursive descendant walker (internal nodes + tips)
# -----------------------------
get_all_descendants <- function(tree, start_node) {
  edge <- tree$edge
  out <- integer(0)
  stack <- start_node

  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]

    children <- edge[edge[, 1] == current, 2]
    if (length(children) == 0) next

    out <- c(out, children)
    stack <- c(children, stack)
  }

  unique(out)
}

# -----------------------------
# 5) Build binary node_scores for internal nodes AND tip labels
# -----------------------------
n_tip <- ape::Ntip(tree_trim2)

internal_scores <- setNames(
  rep(0L, tree_trim2$Nnode),
  as.character((n_tip + 1):(n_tip + tree_trim2$Nnode))
)

tip_scores <- setNames(
  rep(0L, n_tip),
  tree_trim2$tip.label
)

# Paint larger clades first so parent clades don't get overwritten by smaller ones
tmp$subtree_size <- vapply(tmp$clade_node, function(nd) {
  length(get_all_descendants(tree_trim2, nd))
}, numeric(1))

tmp <- tmp[order(tmp$subtree_size, decreasing = TRUE), ]

for (i in seq_len(nrow(tmp))) {
  mrca_node <- tmp$clade_node[i]

  desc <- get_all_descendants(tree_trim2, mrca_node)
  desc_internal <- desc[desc > n_tip]
  desc_tips <- desc[desc <= n_tip]

  internal_scores[as.character(c(mrca_node, desc_internal))] <- 1L
  if (length(desc_tips) > 0) {
    tip_scores[tree_trim2$tip.label[desc_tips]] <- 1L
  }
}

node_scores <- c(internal_scores, tip_scores)

cat("Enriched scored:", sum(node_scores == 1L, na.rm = TRUE), "\n")
cat("Unenriched scored:", sum(node_scores == 0L, na.rm = TRUE), "\n")
cat("Total scored entries:", length(node_scores), "\n")

# -----------------------------
# 6) Plot
# -----------------------------
pdf(
  "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/inflam_tree_HipnBlood_Surgical_binary.pdf",
  width = 10, height = 8
)

plot_phylo_heatmap2(
  tree_trim2,
  df_var = NULL,
  cell_annot = list(clone_df, cell_type_df2),
  annot_pal = list(
    clone_pal,
    c(
      P = "white",
      H = "white",
      NK = "white",
      Mono = "#ffffff",
      B_cell = "white",
      T_cell = "white",
      Oligos = "white",
      HOMEO = "#ffffff",
      INFLAM = "#c410db",
      `T cell` = "white"
    )
  ),
  node_scores = node_scores,
  node_score_limits = c(0, 1),
  layered = FALSE,
  node_conf = FALSE,
  dot_size = 2,
  branch_width = 0.35,
  branch_length = FALSE,
  show_variant_names = FALSE,
  title = "Positively enriched clades on trimmed tree"
)



#####       PATH analysis     ######
##### PATH analysis ######
##### PATH analysis ######
library(PATH)
library(ape)
library(Matrix)
library(ggplot2)
library(phytools)

ord <- c("B cell", "CD4", "CD8", "Mystery?", "Microglia", "Monocytes", "Oligo", "P Mcg/PVM?", "T cell")

tree_path <- tree_trim
tree_path$edge.length <- rep(1, nrow(tree_path$edge))
tree_path <- phytools::midpoint.root(tree_path)
tree_path <- ape::multi2di(tree_path)
tree_path$edge.length <- rep(1, nrow(tree_path$edge))

cell_type_df$cell <- as.character(cell_type_df$cell)
cell_type_df$annot <- as.character(cell_type_df$annot)
cell_type_df$annot[cell_type_df$annot == "" | is.na(cell_type_df$annot)] <- NA
cell_type_df <- cell_type_df[!is.na(cell_type_df$annot), ]

cat("\n--- checkpoint 1: raw annotations ---\n")
print(table(cell_type_df$annot, useNA = "ifany"))
cat("unique annot:\n")
print(sort(unique(cell_type_df$annot)))

cat("\nord used:\n")
print(ord)

common <- intersect(tree_path$tip.label, cell_type_df$cell)
cat("\n--- checkpoint 2: overlap ---\n")
cat("tree tips before prune:", length(tree_path$tip.label), "\n")
cat("df rows before prune:", nrow(cell_type_df), "\n")
cat("common cells:", length(common), "\n")

tree_path <- drop.tip(tree_path, setdiff(tree_path$tip.label, common))
cell_type_df <- cell_type_df[cell_type_df$cell %in% common, ]
cell_type_df <- cell_type_df[match(tree_path$tip.label, cell_type_df$cell), ]

cat("\n--- checkpoint 3: alignment after prune/match ---\n")
cat("tree tips:", length(tree_path$tip.label), "\n")
cat("df rows:", nrow(cell_type_df), "\n")
cat("same order:", all(tree_path$tip.label == cell_type_df$cell), "\n")
cat("any duplicated cells in df:", anyDuplicated(cell_type_df$cell), "\n")

cell_type_df$annot <- factor(cell_type_df$annot, levels = ord)

cat("\n--- checkpoint 4: factor counts ---\n")
print(table(cell_type_df$annot, useNA = "ifany"))
cat("levels(cell_type_df$annot):\n")
print(levels(cell_type_df$annot))

W <- ape::vcv(tree_path)
diag(W) <- 0
rs <- rowSums(W)

cat("\n--- checkpoint 5: W before normalization ---\n")
cat("dim(W):", paste(dim(W), collapse = " x "), "\n")
cat("range(rowSums(W)) before norm:", paste(range(rs), collapse = " to "), "\n")
cat("any zero row sums:", any(rs == 0), "\n")

W <- W / rs
W <- as.matrix(W)

cat("\n--- checkpoint 6: W after normalization ---\n")
cat("dim(W):", paste(dim(W), collapse = " x "), "\n")
cat("range(rowSums(W)) after norm:", paste(range(rowSums(W)), collapse = " to "), "\n")
cat("rownames match tree tips:", all(rownames(W) == tree_path$tip.label), "\n")

X <- model.matrix(~ 0 + annot, data = cell_type_df)
colnames(X) <- sub("^annot", "", colnames(X))
X <- X[, ord, drop = FALSE]
X <- as.matrix(X)
rownames(X) <- cell_type_df$cell

cat("\n--- checkpoint 7: X matrix ---\n")
cat("dim(X):", paste(dim(X), collapse = " x "), "\n")
cat("colSums(X):\n")
print(colSums(X))
cat("row sums of X, unique values:\n")
print(unique(rowSums(X)))
cat("rownames(X) match rownames(W):", all(rownames(X) == rownames(W)), "\n")

stopifnot(nrow(X) == nrow(W))
stopifnot(all(rownames(X) == rownames(W)))

out <- xcor(X, W)
zmat <- out$Z.score
rownames(zmat) <- ord
colnames(zmat) <- ord

cat("\n--- checkpoint 8: zmat ---\n")
cat("dim(zmat):", paste(dim(zmat), collapse = " x "), "\n")
cat("all.equal(zmat, t(zmat)):\n")
print(all.equal(zmat, t(zmat)))
cat("max abs asymmetry:\n")
print(max(abs(zmat - t(zmat))))
cat("diagonal:\n")
print(diag(zmat))
cat("range(zmat):", paste(range(zmat), collapse = " to "), "\n")

zmat_plot <- (zmat + t(zmat)) / 2

cat("\n--- checkpoint 9: symmetrized zmat ---\n")
cat("all.equal(zmat_plot, t(zmat_plot)):\n")
print(all.equal(zmat_plot, t(zmat_plot)))
cat("rounded symmetrized zmat:\n")
print(round(zmat_plot, 2))

plotdat <- expand.grid(i = seq_along(ord), j = seq_along(ord))
plotdat$xlab <- ord[plotdat$i]
plotdat$ylab <- ord[plotdat$j]
plotdat$z_score <- mapply(function(a, b) zmat_plot[a, b], plotdat$xlab, plotdat$ylab)
plotdat$z_score[plotdat$j > plotdat$i] <- NA

plotdat$x <- factor(plotdat$xlab, levels = ord)
plotdat$y <- factor(plotdat$ylab, levels = rev(ord))

cat("\n--- checkpoint 10: plotting data ---\n")
cat("n plot rows:", nrow(plotdat), "\n")
cat("non-NA tiles:", sum(!is.na(plotdat$z_score)), "\n")
cat("expected lower triangle tiles:", length(ord) * (length(ord) + 1) / 2, "\n")

p <- ggplot(plotdat, aes(x = x, y = y, fill = z_score)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_point(data = subset(plotdat, !is.na(z_score)), size = 0.8, color = "black") +
  scale_fill_gradient2(
    low = "#6b5bbd",
    mid = "white",
    high = "#f03b20",
    midpoint = 0,
    na.value = "white"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRunAD_copy/path_heatmap_triangle.pdf", width = 8, height = 7)
print(p)
dev.off()




#surgical PATH analysis
we propagated z-scores from clades to their descendants by assigning each branch the maximum absolute z-score encountered along the path from the root to that node, allowing branches to be colored by the enrichment of their most significant ancestral clade.


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
  cutoff = 0.17
)
dev.off()
# 5. Trim tree based on the confidence threshold of previous plot

tree_trim <- trim_tree(md$tree_annot, conf = 0.145)

pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/Trimmed_tree_Ctx_Surgical.pdf", width = 10, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  node_conf = TRUE,
  dot_size = 2,
  branch_length = FALSE,
  title = "Trimmed tree Blood+Brain_Surgical"
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
pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/Clone_tree_HipnBlood_Surgical.pdf", width = 10, height = 8)
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
cell_type_df$annot <- as.character(cell_type_df$annot)

# Only convert truly empty/NA to "other", but keep Mystery?
cell_type_df$annot[cell_type_df$annot == "" | is.na(cell_type_df$annot)] <- "Unannotated"

# Filter out ONLY the "Unannotated" cells (not Mystery?)
cell_type_df <- cell_type_df[cell_type_df$annot != "Unannotated", ]

# Verify Mystery? is still there
print(table(cell_type_df$annot))
print(sum(cell_type_df$annot == "Mystery?"))
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

#Create tissue type annotation dataframe
df_with_annot <- df_with_annot %>%
  mutate(
    annot = substr(gsub("^(.)(_.*)", "\\1", cell), 1, 1),
    annot = if_else(annot %in% c("P", "H"), annot, "Unknown")  # safety if you ever see non‑P/H
  ) 

#Create cell type annotation dataframe
cell_type_df <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/tissue_mitodrift_df_combined.csv", stringsAsFactors = FALSE)
cell_type_df




                                                    ### clade enrichment analysis ### 
shared <- intersect(tree_trim$tip.label, cell_type_df$cell)
tree_trim2 <- ape::keep.tip(tree_trim, shared)
cell_type_df2 <- cell_type_df[cell_type_df$cell %in% shared, ]
trait <- setNames(cell_type_df2$annot, cell_type_df2$cell)
trait <- trait[tree_trim2$tip.label]
mono <- as.numeric(trait == "INFLAM")
names(mono) <- names(trait)
# binary monocyte trait: 1 = Mono, 0 = not Mono
mono <- as.numeric(trait == "INFLAM")
names(mono) <- tree_trim2$tip.label
# global mean and SD across all tips
global_mean <- mean(mono, na.rm = TRUE)
global_sd <- sd(mono, na.rm = TRUE)
# get internal nodes
Ntip <- length(tree_trim2$tip.label)
nodes <- (Ntip + 1):(Ntip + tree_trim2$Nnode)
# compute clade mean and deviation for each internal node
mono_clade_stats <- data.frame(
  node = nodes,
  size = NA_integer_,
  clade_mean = NA_real_,
  deviation = NA_real_
)
for (i in seq_along(nodes)) {
  node <- nodes[i]
  desc <- phangorn::Descendants(tree_trim2, node, type = "tips")[[1]]
  tip_names <- tree_trim2$tip.label[desc]
  vals <- mono[tip_names]
  
  mono_clade_stats$size[i] <- length(vals)
  mono_clade_stats$clade_mean[i] <- mean(vals, na.rm = TRUE)
  mono_clade_stats$deviation[i] <- mono_clade_stats$clade_mean[i] - global_mean
}
head(mono_clade_stats)

#Deviation from gloabl mean / size-adjusted SE = z score, use clone_df to get clade size
mono_clade_stats <- clone_df %>%
  group_by(clade) %>%
  summarise(
    size = n(),
    clade_mean = mean(mono[cell], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    deviation = clade_mean - global_mean,
    se = global_sd / sqrt(size),
    z = deviation / se
  )
head(mono_clade_stats)

#calculate p value from z score 
min_clade_size <- 10 #select major clades, do this per cell type
mono_clade_stats$p <- NA_real_
keep <- !is.na(mono_clade_stats$z) & mono_clade_stats$size >= min_clade_size
# two-sided p-values from standard normal
mono_clade_stats$p[keep] <- 2 * pnorm(-abs(mono_clade_stats$z[keep]))
# Benjamini-Hochberg correction
mono_clade_stats$padj <- NA_real_
mono_clade_stats$padj[keep] <- p.adjust(mono_clade_stats$p[keep], method = "BH")
head(mono_clade_stats)


#'clade enrichment'
library(ape)
library(dplyr)

# -----------------------------
# 1) Prepare trimmed tree and matching cell annotation
# -----------------------------
shared <- intersect(tree_trim$tip.label, cell_type_df$cell)
tree_trim2 <- ape::keep.tip(tree_trim, shared)
cell_type_df2 <- cell_type_df[cell_type_df$cell %in% shared, ]

# -----------------------------
# 2) Use existing mono_clade_stats, keep only positive enrichment
#    Expected columns: clade, z, size, padj
# -----------------------------
mono_clade_stats2 <- mono_clade_stats %>%
  distinct(clade, .keep_all = TRUE) %>%
  mutate(clade = as.character(clade)) %>%
  filter(size >= 10, z > 0)
# Optional stricter version:
# mono_clade_stats2 <- mono_clade_stats %>%
#   distinct(clade, .keep_all = TRUE) %>%
#   mutate(clade = as.character(clade)) %>%
#   filter(size >= 10, z > 0, !is.na(padj), padj < 0.05)

# -----------------------------
# 3) Map each clade to its MRCA node in the trimmed tree
# -----------------------------
clade_tips <- clone_df %>%
  group_by(clade) %>%
  summarise(
    tips = list(intersect(unique(cell), tree_trim2$tip.label)),
    .groups = "drop"
  )

clade_node_map2 <- lapply(seq_len(nrow(clade_tips)), function(i) {
  tips_i <- clade_tips$tips[[i]]
  if (length(tips_i) < 2) return(NULL)

  node_i <- ape::getMRCA(tree_trim2, tips_i)
  if (is.null(node_i) || is.na(node_i)) return(NULL)

  data.frame(
    clade = as.character(clade_tips$clade[i]),
    clade_node = as.integer(node_i),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

tmp <- merge(
  clade_node_map2,
  mono_clade_stats2[, c("clade", "z", "size", "padj")],
  by = "clade",
  all.x = TRUE
)

tmp <- tmp[!is.na(tmp$z), ]
tmp <- tmp[tmp$size >= 10, ]

cat("Positive enriched clades mapped:", nrow(tmp), "\n")

# -----------------------------
# 4) Recursive descendant walker (internal nodes + tips)
# -----------------------------
get_all_descendants <- function(tree, start_node) {
  edge <- tree$edge
  out <- integer(0)
  stack <- start_node

  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]

    children <- edge[edge[, 1] == current, 2]
    if (length(children) == 0) next

    out <- c(out, children)
    stack <- c(children, stack)
  }

  unique(out)
}

# -----------------------------
# 5) Build binary node_scores for internal nodes AND tip labels
# -----------------------------
n_tip <- ape::Ntip(tree_trim2)

internal_scores <- setNames(
  rep(0L, tree_trim2$Nnode),
  as.character((n_tip + 1):(n_tip + tree_trim2$Nnode))
)

tip_scores <- setNames(
  rep(0L, n_tip),
  tree_trim2$tip.label
)

# Paint larger clades first so parent clades don't get overwritten by smaller ones
tmp$subtree_size <- vapply(tmp$clade_node, function(nd) {
  length(get_all_descendants(tree_trim2, nd))
}, numeric(1))

tmp <- tmp[order(tmp$subtree_size, decreasing = TRUE), ]

for (i in seq_len(nrow(tmp))) {
  mrca_node <- tmp$clade_node[i]

  desc <- get_all_descendants(tree_trim2, mrca_node)
  desc_internal <- desc[desc > n_tip]
  desc_tips <- desc[desc <= n_tip]

  internal_scores[as.character(c(mrca_node, desc_internal))] <- 1L
  if (length(desc_tips) > 0) {
    tip_scores[tree_trim2$tip.label[desc_tips]] <- 1L
  }
}

node_scores <- c(internal_scores, tip_scores)

cat("Enriched scored:", sum(node_scores == 1L, na.rm = TRUE), "\n")
cat("Unenriched scored:", sum(node_scores == 0L, na.rm = TRUE), "\n")
cat("Total scored entries:", length(node_scores), "\n")

# -----------------------------
# 6) Plot
# -----------------------------
pdf(
  "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/inflam_tree_HipnBlood_Surgical_binary.pdf",
  width = 10, height = 8
)

plot_phylo_heatmap2(
  tree_trim2,
  df_var = NULL,
  cell_annot = list(clone_df, cell_type_df2),
  annot_pal = list(
    clone_pal,
    c(
      P = "white",
      H = "white",
      NK = "white",
      Mono = "#ffffff",
      B_cell = "white",
      T_cell = "white",
      Oligos = "white",
      HOMEO = "#ffffff",
      INFLAM = "#c410db",
      `T cell` = "white"
    )
  ),
  node_scores = node_scores,
  node_score_limits = c(0, 1),
  layered = FALSE,
  node_conf = FALSE,
  dot_size = 2,
  branch_width = 0.35,
  branch_length = FALSE,
  show_variant_names = FALSE,
  title = "Positively enriched clades on trimmed tree"
)



#####       PATH analysis     ######
##### PATH analysis 220626#
# keep only P_1235 and C_5
mut2 <- mut_dat[mut_dat$annot %in% c("P_1235", "C_5"), ]
mut2$group <- factor(mut2$annot, levels = c("P_1235", "C_5"))

# align to tree
tree2 <- tree_trim
tree2$edge.length <- rep(1, nrow(tree2$edge))
tree2 <- phytools::midpoint.root(tree2)
tree2 <- ape::multi2di(tree2)
tree2$edge.length <- rep(1, nrow(tree2$edge))

common <- intersect(tree2$tip.label, mut2$cell)
tree2 <- drop.tip(tree2, setdiff(tree2$tip.label, common))
mut2 <- mut2[mut2$cell %in% common, ]
mut2 <- mut2[match(tree2$tip.label, mut2$cell), ]

# W matrix
W <- ape::vcv(tree2)
diag(W) <- 0
W <- W / rowSums(W)
W <- as.matrix(W)

# X matrix
X <- model.matrix(~ 0 + group, data = mut2)
colnames(X) <- sub("^group", "", colnames(X))
X <- as.matrix(X)
X <- X[, c("P_1235", "C_5"), drop = FALSE]
rownames(X) <- mut2$cell

# PATH test
out <- xcor(X, W)

zmat <- out$Z.score
print(zmat)
print(out$phy_cor)
print(out$one.sided.pvalue)