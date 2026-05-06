library(mitodrift)
library(data.table)

base_dir    <- "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/MonoSubsets"
results_dir <- file.path(base_dir, "Results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

donor_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
n_written  <- 0

for (donor_dir in donor_dirs) {
  group_name <- basename(donor_dir)

  sub_dirs <- list.dirs(donor_dir, recursive = FALSE, full.names = TRUE)
  sub_dirs <- sub_dirs[grepl("^subsample_seed[0-9]+$", basename(sub_dirs))]

  for (sub_dir in sub_dirs) {
    subsample_name <- basename(sub_dir)

    mut_file <- file.path(sub_dir, "mut_dat_sub.csv")
    md_file  <- file.path(sub_dir, "mitodrift_object.rds")

    # Skip missing subsamples
    if (!file.exists(mut_file) || !file.exists(md_file)) {
      warning("Missing mut_dat_sub.csv or mitodrift_object.rds in: ", sub_dir)
      next
    }

    md        <- readRDS(md_file)
    tree_trim <- trim_tree(md$tree_annot, conf = 0.01138)
    clone_df  <- assign_clones_polytomy(tree_trim)

    # Safety check
    if (!is.data.frame(clone_df) ||
        !all(c("cell", "clade", "clade_node", "annot", "size", "frac") %in% names(clone_df))) {
      warning("clone_df has wrong structure in: ", sub_dir)
      next
    }

    # Sanity check: warn if one big clone dominates
    top_frac <- max(clone_df$frac)
    if (top_frac > 0.8) {
      warning("One big clone >80% in ", sub_dir, " (top_frac = ", round(top_frac, 3), ")")
    }

    # Save to /MonoSubsets/Results as requested
    out_file <- file.path(
      results_dir,
      paste0(group_name, "_", subsample_name, "_clone_df.csv")
    )
    fwrite(clone_df, out_file)
    n_written <- n_written + 1
  }
}

message("Wrote ", n_written, " clone_df CSVs")


library(ggplot2)
library(dplyr)

# Assuming plot_data is already calculated as in your snippet
p <- ggplot(plot_data, aes(x = donor, y = mtSDI_mean)) +
  geom_point(size = 2.5, color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.1, linewidth = 0.3, color = "black") +
  theme_bw(base_size = 12) +
  theme(
    # Disable vertical grid lines to remove the "separated" look
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "mtSDI") +
  # Shrink the white space around the points to bring them closer
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2)))

print(p)
ggsave("mtSDI_dotplot_compact.pdf", p, width = 1.5, height = 3.0)