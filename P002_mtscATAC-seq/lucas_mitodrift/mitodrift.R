library(mitodrift)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtree)
library(ape)

# 1. Your mutation data (same one you used as input)
mut_dat <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex0.05/Cortex_surgical_long_variant_matrix.csv")

# 2. Your model output
md <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex0.05/mitodrift_object.rds")

tree <- md$tree_annot

#cehck tree structure 
names(md)

mcmc_trace <- qread("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex/tree_ml_trace.qs2")
tail(mcmc_trace$asdsf)  # Last ASDSF values
diag <- readRDS("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex/tree_mcmc_diag.rds")
diag$asdsf  # Full ASDSF trajectory




# 3. Trim the tree
tree_trim <- trim_tree(tree, conf = 0.2)

# 4. Plot
pdf("mitodriftHippocamp_plot.pdf", width = 20, height = 16)
plot_phylo_heatmap2(
  tree,
  mut_dat,
  dot_size = 1,
  branch_width = 0.3,
  branch_length = TRUE,
  node_conf = TRUE,
  ladderize = TRUE,
  ylim = NULL,
  het_max = 1,
  node_lab = FALSE,
  title = "Hippocamp_plot"
)
dev.off()



tree_trim <- trim_tree(tree, conf = 0.1, collapse_trivial = TRUE)



pdf("Trimmed_tree_0.1_PMBrain.pdf", width = 15, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  node_conf = FALSE,
  dot_size = 2,
  branch_length = TRUE,
  conf_label = FALSE,
  tip_lab = FALSE,
  title = "Trimmed tre"
)
dev.off()
# 6. Save the variants in the trimmed tree
tree_variants = unique(mut_dat$variant)
write.csv(tree_variants, "tree_variants_PMbrain.csv", row.names = FALSE)
#head(mut_dat, n=20)
# 7. Assign clones


clone_df <- assign_clones_polytomy(tree_trim)


# Remove clades where size == 1 (only one cell in the clade)
clone_df_nosingletons <- clone_df %>%
  filter(size > 1)
# Check how many cells and clades you kept
cat("Number of cells in clades with size > 1:", nrow(clone_df_nosingletons), "\n")
cat("Number of unique clades with size > 1:", n_distinct(clone_df_nosingletons$clade), "\n")


# 1) Get singleton cells (cells in clades with size == 1)
singleton_cells <- clone_df %>%
  filter(size == 1) %>%
  pull(cell)
cat("Number of singleton cells (size == 1):", length(singleton_cells), "\n")
# Prune them from the tree
tree_trim_nosingletons <- drop.tip(
  phy = tree_trim,
  tip = intersect(singleton_cells, tree_trim$tip.label),
  collapse.singles = TRUE
)
mut_dat_nosingletons <- mut_dat %>%
  filter(!cell %in% singleton_cells)
  # Sanity check: tree tips match mut_dat cells
stopifnot(
  setequal(
    sort(tree_trim_nosingletons$tip.label),
    sort(unique(mut_dat_nosingletons$cell))
  )
)
clone_df_nosingletons <- clone_df %>%
  filter(size > 1)

# 8. Visualise the clones
clade_order <- unique(clone_df$clade)
clone_pal <- make_clade_pal(length(clade_order), labels = clade_order,
                            pal = "Dark2", cycle_len = 8, cycle_shift = 0)
pdf("Clone_tree_run2_surgicalBrain.pdf", width = 15, height = 8)
plot_phylo_heatmap2(
  tree_trim,
  mut_dat,
  cell_annot = clone_df,
  annot_pal = clone_pal,
  node_conf = FALSE,
  conf_label = FALSE,
  dot_size = 2,
  branch_length = TRUE,
  mark_low_cov = FALSE,
  node_lab = FALSE,
  title = "Clones on trimmed tree"
)
dev.off()
head(clone_df)
head(clone_pal)
pdf("Circ_chart_run2_PMBrain.pdf", width = 10, height = 8)
plot_phylo_circ(
  tree_trim_nosingletons,
  cell_annot = clone_df_nosingletons,
  annot_pal = clone_pal,
  annot_legend = FALSE,
  title = "Circular layout"
)
dev.off()


library(dplyr)
combined_annot <- clone_df_nosingletons %>%
  left_join(leiden_simple, by = "cell") %>%
  select(cell, clade, leiden)



# Palette for BOTH bars
leiden_levels <- sort(unique(celltypes_df$leiden))
n_leiden <- length(leiden_levels)
leiden_colors <- RColorBrewer::brewer.pal(min(12, n_leiden), "Set3")
if (n_leiden > 12) leiden_colors <- colorRampPalette(leiden_colors)(n_leiden)

# Name them by leiden numbers
names(leiden_colors) <- leiden_levels

# Combined palette
clone_pal_extended <- list(
  clade = clone_pal,      # Your existing clone colors
  leiden = leiden_colors  # Auto-generated leiden colors
)
celltypes_df <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex0.05/clades_with_leiden.csv")

library(ggplot2)

pr_df <- compute_variant_pr_curve(tree, mut_dat)

p1 <- plot_prec_recall_vs_conf(pr_df, 
                               sample_name = "Surgical_Cortex Variant PR",
                               cutoff = 0.15)

#Save
ggsave("variant_PR_curve_Surgical_Cortex.pdf", p1, width=8, height=6)




tree$node.label
summary(as.numeric(tree$node.label))  # Min/Max/Quartiles
hist(as.numeric(tree$node.label), breaks=20, main="Node Confidences")  # Distribution
table(as.numeric(tree$node.label) < 0.125)  # Count below common cutoff


pr_df


 write.csv(clone_df, "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex0.05/clone_assignments_wsingletons.csv", row.names = FALSE)


library(dplyr)


library(data.table)

# Paths
cortex_file <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Cat0.05/Cortex_surgical_long_variant_matrix.csv"
blood_file <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Cat0.05/Fullblood_surgical_long_variant_matrix.csv"
output_path <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/Cat0.05/combined_blood_cortex_long_variant_matrix.csv"

# Load CSV (comma-separated, no headers)
cat("Loading CSV files...\n")
cortex_dt <- fread(cortex_file, sep=",", header=FALSE, col.names=c("cell_id", "variant_id", "cortex_count", "blood_count"))
blood_dt <- fread(blood_file, sep=",", header=FALSE, col.names=c("cell_id", "variant_id", "cortex_count", "blood_count"))

cat("Shapes:", dim(cortex_dt), dim(blood_dt), "\n")

# Prefix C_ for cortex, B_ for blood
cat("Prefixing...\n")
cortex_dt[, cell_id := paste0("C_", cell_id)]
cortex_dt[, variant_id := paste0("C_", variant_id)]
blood_dt[, cell_id := paste0("B_", cell_id)]
blood_dt[, variant_id := paste0("B_", variant_id)]

# Combine
combined_dt <- rbind(cortex_dt, blood_dt)

# Inspect
cat("Combined rows:", nrow(combined_dt), "\n")
print(head(combined_dt, 6))

# Save EXACT SAME FORMAT: CSV, comma, NO headers
fwrite(combined_dt, output_path, sep=",", col.names=FALSE, row.names=FALSE)
cat("Saved identical format:", output_path, "\n")
cat("Size:", round(file.size(output_path)/1e9, 2), "GB\n")


leiden_df <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgicalCortex0.05/clades_with_leiden.csv")
leiden_simple <- leiden_df[, c(1, 7)]
leiden_simple

library(dplyr)
library(RColorBrewer)

library(dplyr)
library(RColorBrewer)

library(dplyr)
library(RColorBrewer)

## 1. Join leiden_simple into clone_df_nosingletons
clone_df_nosingletons <- clone_df_nosingletons %>%
  left_join(
    leiden_simple,
    by = "cell"
  )

## 2. Create leiden_simple in clone_df_nosingletons
if ("leiden" %in% colnames(clone_df_nosingletons)) {
  clone_df_nosingletons$leiden_simple <- clone_df_nosingletons$leiden
}

## 3. Sanity check
if (!"leiden_simple" %in% colnames(clone_df_nosingletons)) {
  stop("leiden_simple not found in clone_df_nosingletons; check join key.")
}
if (all(is.na(clone_df_nosingletons$leiden_simple))) {
  stop("All leiden_simple values are NA; check that leiden_simple cells match tree tips.")
}

## 4. Prepare clade annotation
clade_order <- unique(clone_df_nosingletons$clade)

clone_df_nosingletons$clade <- factor(
  clone_df_nosingletons$clade,
  levels = clade_order
)

clone_pal <- make_clade_pal(
  length(clade_order),
  labels = clade_order,
  pal = "Dark2",
  cycle_len = 8,
  cycle_shift = 0
)

## 5. Prepare leiden_simple annotation
leiden_simple_order <- sort(unique(clone_df_nosingletons$leiden_simple[!is.na(clone_df_nosingletons$leiden_simple)]))

clone_df_nosingletons$leiden_simple <- factor(
  clone_df_nosingletons$leiden_simple,
  levels = leiden_simple_order
)

if (length(leiden_simple_order) <= 8) {
  leiden_cols <- brewer.pal(8, "Set2")[seq_along(leiden_simple_order)]
} else {
  leiden_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(leiden_simple_order))
}

leiden_pal <- setNames(leiden_cols, leiden_simple_order)

## 6. Build wide annotation data frame (NOT a list)
cell_annot_plot <- clone_df_nosingletons %>%
  dplyr::select(
    cell,
    clade,
    leiden_simple
  ) %>%
  arrange(cell)

## 7. Combine palettes as named list keyed by column names
annot_pal <- list(
  clade         = clone_pal,
  leiden_simple = leiden_pal
)

## 8. Plot with layered = TRUE
pdf("Clone_tree_run2_PMBrain.pdf", width = 15, height = 8)

p <- plot_phylo_heatmap2(
  phylo            = tree_trim_nosingletons,
  df_var           = mut_dat_nosingletons,
  cell_annot       = cell_annot_plot,
  annot_pal        = annot_pal,
  layered          = TRUE,
  node_conf        = FALSE,
  conf_label       = FALSE,
  dot_size         = 2,
  branch_length    = TRUE,
  mark_low_cov     = FALSE,
  node_lab         = FALSE,
  annot_legend     = TRUE,
  annot_bar_height = 0.1,
  title            = "Clades and Leiden clusters on trimmed tree"
)

print(p)

dev.off()