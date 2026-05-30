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
library(PATH)
library(ape)
library(Matrix)
library(ggplot2)
library(phytools)

ord <- c("B cell", "NK", "CD8 T cell", "CD4 T cell", "Mono", "T cell", "INFLAM", "HOMEO", "Oligos")

tree_path <- tree_trim
tree_path$edge.length <- rep(1, nrow(tree_path$edge))
tree_path <- phytools::midpoint.root(tree_path)
tree_path <- ape::multi2di(tree_path)
tree_path$edge.length <- rep(1, nrow(tree_path$edge))

cell_type_df$cell <- as.character(cell_type_df$cell)
cell_type_df$annot <- as.character(cell_type_df$annot)
cell_type_df$annot[cell_type_df$annot == "" | is.na(cell_type_df$annot)] <- "other"
cell_type_df <- cell_type_df[cell_type_df$annot != "other", ]

common <- intersect(tree_path$tip.label, cell_type_df$cell)
tree_path <- drop.tip(tree_path, setdiff(tree_path$tip.label, common))
cell_type_df <- cell_type_df[cell_type_df$cell %in% common, ]
cell_type_df <- cell_type_df[match(tree_path$tip.label, cell_type_df$cell), ]

cell_type_df$annot <- factor(cell_type_df$annot, levels = ord)

W <- ape::vcv(tree_path)
diag(W) <- 0
W <- W / rowSums(W)

X <- model.matrix(~ 0 + annot, data = cell_type_df)
colnames(X) <- sub("^annot", "", colnames(X))
X <- X[, ord, drop = FALSE]

out <- xcor(X, W)
zmat <- out$Z.score
rownames(zmat) <- ord
colnames(zmat) <- ord

plotdat <- expand.grid(i = seq_along(ord), j = seq_along(ord))
plotdat$xlab <- ord[plotdat$i]
plotdat$ylab <- ord[plotdat$j]
plotdat$z_score <- mapply(function(a, b) zmat[a, b], plotdat$xlab, plotdat$ylab)

plotdat <- plotdat[plotdat$j <= plotdat$i, ]

plotdat$x <- plotdat$i
plotdat$y <- length(ord) - plotdat$j + 1

p <- ggplot(plotdat, aes(x = x, y = y, fill = z_score)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_point(size = 0.8, color = "black") +
  scale_fill_gradient2(low = "#6b5bbd", mid = "white", high = "#f03b20", midpoint = 0) +
  scale_x_continuous(
    breaks = seq_along(ord),
    labels = ord,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(ord),
    labels = rev(ord),
    expand = c(0, 0)
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

pdf("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/CatRun210526_copy/path_heatmap_triangle.pdf", width = 8, height = 7)
print(p)
dev.off()