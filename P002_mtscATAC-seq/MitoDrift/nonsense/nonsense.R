suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(stringr)
    library(dplyr)
    library(parallel)
    library(glue)
    library(patchwork)
    library(ape)
    library(igraph)
    library(ggtree)
    library(ggraph)
    library(ggpubr)
    library(tidygraph)
    library(ggrepel)
    library(ggtreeExtra)
    library(mitodrift)
    library(ggrastr)
})

filter = dplyr::filter
options(repr.matrix.max.cols=15, repr.matrix.max.rows=50)
rename = dplyr::rename
home='/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/HPC_outputs/PM_copy'
R.utils::sourceDirectory(file.path(home, "R"))

sample_id = 'PM_cat'
phy_annot = read.tree(glue('/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/HPC_outputs/PM_copy/tree_annotated.newick'))
mut_dat = fread(glue('/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/HPC_outputs/PM_copy/PM_combined_matrix.csv'))
cell_annot <- read.csv("/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/nonsense/cell_with_prefix_corrected.csv")

options(repr.plot.width = 2.75, repr.plot.height = 2, repr.plot.res = 250)

library(mitodrift)

p_conf = 0.1
phy_trim = phy_annot %>% trim_tree(conf = p_conf)

pr_df_var = compute_variant_pr_curve(phy_annot, mut_dat, j_thres = 0.5, min_vaf = 0.05, ncores = 8)
p = plot_prec_recall_vs_conf(pr_df_var, cutoff = p_conf, legend = F) + 
    ggtitle(glue('PM_cat - Variants'))

p

ggsave(glue('/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/nonsense/variant_diag_curve.pdf'), p, width = 2.75, height = 2)




library(mitodrift)
options(repr.plot.width = 6.75, repr.plot.height = 4, repr.plot.res = 250)
source("/mnt/claw-raid/elliot/P002_mtscATAC-seq/scripts/path_state_fate_analysis/R/utils.R")

pal_ct <- get_discrete_colors(cell_annot$prefix, "Set3")

p = plot_phylo_heatmap2(
    phy_trim,
    mut_dat %>% group_by(variant) %>% filter(sum(a>0)>=10),
    node_conf = F,
    branch_length = F,
    branch_width = 0.05,
    het_max = 0.2,
    annot_legend = T,
    cell_annot = list(
        'Celltype' = cell_annot %>% mutate(annot = prefix)
    ),
    annot_pal = list(pal_ct),
    annot_title_size = 0,
    show_variant_names = FALSE
) %>% suppressWarnings() %>%
suppressMessages()


ggsave(
    "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/nonsense/tree.pdf",
    p,
    width = 6.75,
    height = 4
)


phy_trim <- drop_singletons(phy_trim)


options(repr.plot.width = 3, repr.plot.height = 2, repr.plot.res = 200)

cluster_dict = cell_annot %>% {setNames(.$prefix, .$cell)}

source("/mnt/claw-raid/elliot/P002_mtscATAC-seq/scripts/path_state_fate_analysis/R/lineage_coupling.R")

state_dict <- cluster_dict[phy_trim$tip.label]

res_path <- run_path (phy_trim, state_dict = state_dict, model = "sibling", norm = TRUE, min_count = 10, diag_only = FALSE)



p = plot_coupling_triangle(res_path, mark_signif = TRUE, limits = c(-2.5, 2.5), statistic = "z",
    signif_only = FALSE, q_thres = 0.2, rev_order = FALSE, highlight = NULL, highlight_color = "firebrick")
 +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right",
    ) )

ggsave(
    "/mnt/claw-raid/elliot/P002_mtscATAC-seq/MitoDrift/nonsense/triangle.pdf",
    p,
    width = 6.75,
    height = 4
)










library(stringr)

# --- ensure consistent state order (IMPORTANT) ---
state_order <- c(
  "P_",
  "C_0","C_1","C_2","C_3","C_4","C_5"
)
# --- node colors: P vs C ---
node_fill_dict <- setNames(
  ifelse(str_detect(state_order, "^P_"), "orange", "brown"),
  state_order
)

# --- node sizes (IMPORTANT FIX) ---
# cluster_dict is cell -> state, so we size by frequency of states

node_size_dict <- table(cluster_dict[names(cluster_dict) %in% names(cluster_dict)])
node_size_dict <- log(as.numeric(node_size_dict) + 1)
names(node_size_dict) <- names(table(cluster_dict))
p <- plot_coupling_panel(
    res_path,
    tmat,
    state_order,
    linewidth_range = c(1, 2),
    node_fill_map = node_fill_dict,
    node_size_map = node_size_dict,
    limits = c(-50, 50),
    self_coupling_limits = c(0, 0.35),
    statistic = "z",
    mark_signif = TRUE,
    panel_spacing = 0.3
) &
guides(
    fill = guide_colourbar(
        title.position = "top",
        title.hjust = 1,
        label.position = "left",
        label.hjust = 1,
        title = "Coupling\nZ-score"
    ),
    colour = guide_colourbar(
        title.position = "top",
        title.hjust = 1,
        label.position = "left",
        label.hjust = 1,
        title = "Transition\nStrength"
    )
)