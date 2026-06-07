library(mitodrift)
library(data.table)
library(ggplot2)

base_dir <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes"
results_dir <- file.path(base_dir, "Results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Subtype folders expected directly under base_dir
subtype_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subtype_dirs <- subtype_dirs[basename(subtype_dirs) %in% c("HOMEO", "INFLAM")]

if (length(subtype_dirs) == 0) {
  stop("No HOMEO/INFLAM directories found under: ", base_dir)
}

mtSDI_results <- data.table()
n_processed <- 0

for (subtype_dir in subtype_dirs) {
  subtype <- basename(subtype_dir)
  cat("Processing subtype:", subtype, "\n")
  
  sub_dirs <- list.dirs(subtype_dir, recursive = FALSE, full.names = TRUE)
  sub_dirs <- sub_dirs[grepl("^subsample_seed[0-9]+$", basename(sub_dirs))]
  
  if (length(sub_dirs) == 0) {
    warning("No subsample_seed dirs found in: ", subtype_dir)
    next
  }
  
  for (sub_dir in sub_dirs) {
    subsample <- basename(sub_dir)
    md_file <- file.path(sub_dir, "mitodrift_object.rds")
    
    if (!file.exists(md_file)) {
      warning("Missing mitodrift_object.rds in: ", sub_dir)
      next
    }
    
    md <- readRDS(md_file)
    tree_trim <- trim_tree(md$tree_annot, conf = 0.125)
    clone_df <- assign_clones_polytomy(tree_trim)
    
    if (!is.data.frame(clone_df)) {
      warning("clone_df is not a data.frame in: ", sub_dir)
      next
    }
    
    clone_df <- as.data.table(clone_df)
    
    req_cols <- c("cell", "clade", "clade_node", "annot", "size", "frac")
    missing_cols <- setdiff(req_cols, names(clone_df))
    if (length(missing_cols) > 0) {
      warning("Missing columns in clone_df for ", sub_dir, ": ", paste(missing_cols, collapse = ", "))
      next
    }
    
    clone_df[, subtype := subtype]
    clone_df[, subsample := subsample]
    
    # Save raw clone assignments for each subsample
    fwrite(
      clone_df,
      file.path(results_dir, paste0(subtype, "_", subsample, "_clone_df.csv"))
    )
    
    # Collapse to unique clones within this subsample
    clone_summary <- unique(clone_df[, .(clade, clade_node, annot, size, frac)])
    
    # Safety: only positive clone fractions
    clone_summary <- clone_summary[frac > 0]
    
    if (nrow(clone_summary) == 0) {
      warning("No positive-fraction clones in: ", sub_dir)
      next
    }
    
    # mtSDI exactly as Shannon diversity on clone fractions
    mtSDI_val <- -sum(clone_summary$frac * log(clone_summary$frac))
    
    mtSDI_results <- rbind(
      mtSDI_results,
      data.table(
        subtype = subtype,
        subsample = subsample,
        mtSDI = mtSDI_val,
        n_clones = nrow(clone_summary),
        top_clone_frac = max(clone_summary$frac),
        total_cells = sum(clone_summary$size)
      ),
      fill = TRUE
    )
    
    n_processed <- n_processed + 1
    cat("  Done:", subsample, "| mtSDI =", round(mtSDI_val, 3), "\n")
  }
}

cat("\nProcessed", n_processed, "subsamples\n")

if (nrow(mtSDI_results) == 0) {
  stop("No mtSDI results generated.")
}

# Summary by subtype
mtSDI_summary <- mtSDI_results[, .(
  mean_mtSDI = mean(mtSDI),
  sd_mtSDI = sd(mtSDI),
  ci_lower = mean(mtSDI) - 1.96 * sd(mtSDI) / sqrt(.N),
  ci_upper = mean(mtSDI) + 1.96 * sd(mtSDI) / sqrt(.N),
  n_subsamples = .N,
  mean_n_clones = mean(n_clones),
  mean_top_clone_frac = mean(top_clone_frac)
), by = subtype]

cat("\nmtSDI summary:\n")
print(mtSDI_summary)

# Optional test between HOMEO and INFLAM
if (length(unique(mtSDI_results$subtype)) == 2) {
  p_val <- wilcox.test(mtSDI ~ subtype, data = mtSDI_results)$p.value
  cat("\nWilcoxon p-value:", format.pval(p_val, digits = 3), "\n")
}

# Save tables
fwrite(mtSDI_results, file.path(results_dir, "mtSDI_by_subsample.csv"))
fwrite(mtSDI_summary, file.path(results_dir, "mtSDI_summary_HOMEO_INFLAM.csv"))

# FIXED plot showing ALL 3 subtypes with proper spacing
p <- ggplot(mtSDI_summary, aes(x = reorder(subtype, -mean_mtSDI), y = mean_mtSDI)) +  # - for descending order
  geom_point(size = 4, color = "black", fill = "white", stroke = 1.2, shape = 21) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.25, linewidth = 0.8, color = "black") +
  geom_text(aes(label = sprintf("n=%d\n%.3f", n_subsamples, mean_mtSDI)), 
            vjust = -1.2, hjust = 0.5, size = 3.5, fontface = "bold", color = "black") +
  coord_cartesian(ylim = c(5.0, 6.1)) +  # Zoom to your data range
  theme_bw(base_size = 14) +
  labs(
    x = "Microglia Subtype",
    y = "mtSDI (mean ± 95% CI)",
    title = "P002: mtSDI by Subtype",
    subtitle = "Monocytes >> Inflammatory > Homeostatic"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40", face = "italic"),
    plot.margin = margin(15, 15, 15, 15)
  )

print(p)
ggsave(file.path(results_dir, "mtSDI_ALL3.pdf"), p, width = 6, height = 5, dpi = 300)



library(data.table)
library(ggplot2)

summary_file <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes/Results/mtSDI_summary_HOMEO_INFLAM.csv"
results_dir   <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/SurgMCGSubtypes/Results"

mtSDI_summary <- fread(summary_file)

# Clean subtype names just in case
mtSDI_summary[, subtype := trimws(subtype)]

print(mtSDI_summary)

p <- ggplot(mtSDI_summary, aes(x = subtype, y = mean_mtSDI)) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 0.4
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "mtSDI",
    title = "mtSDI by subtype"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1)
  )

print(p)

ggsave(
  file.path(results_dir, "mtSDI_summary_from_csv.pdf"),
  p,
  width = 4.5,
  height = 3.5
)