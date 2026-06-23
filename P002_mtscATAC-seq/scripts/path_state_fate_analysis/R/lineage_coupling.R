
##' Run PATH lineage coupling / phylogenetic association test on a tree
##'
##' This wrapper computes pairwise association (Z-score and correlation)
##' between discrete states or continuous features and phylogenetic structure
##' using the PATH package helpers. It accepts either a named vector of
##' categorical states (`state_dict`) keyed by tip names, or a feature matrix
##' (`feature_mat`) with features in rows and cells in columns. The function
##' will subset the tree to the intersection of provided cells and tree tips.
##'
##' Supported models for the tree-distance weighting matrix (`W`):
##' - "bm": Brownian-motion style distance (uses `bm_node_dist`)
##' - "sibling": custom one-node sibling distance (uses `one_node_tree_dist2`)
##' - "inv": inverse tree distance (uses `PATH::inv_tree_dist`)
##'
##' @param tree A rooted `phylo` object.
##' @param state_dict Optional named character vector mapping tip names to discrete states.
##' @param feature_mat Optional numeric matrix of features (rows) × cells (cols). One of `state_dict` or `feature_mat` must be provided.
##' @param model Character; one of "bm", "sibling", or "inv" to choose how to build the distance matrix `W`. Default: "bm".
##' @param norm Logical; whether to normalize the distance matrix. Passed to underlying helpers. Default: TRUE.
##' @param min_count Integer; when `state_dict` is provided, drop states with fewer than `min_count` observations. Default: 0.
##' @param diag_only Logical; if TRUE, return only the diagonal (self) comparisons. Default: FALSE.
##' @return A tibble with columns `state_x`, `state_y`, `z`, `cor`, `p`, and `q` (BH-adjusted p-values). If `diag_only` is TRUE the table will include only the `state` and `z` columns.
##' @examples
##' # Using a small tree and discrete states
##' # res <- run_path(tree, state_dict = my_states, model = "bm")
##'
##' @export
run_path = function(tree, state_dict = NULL, feature_mat = NULL, model = "bm", norm = TRUE, min_count = 0, diag_only = FALSE) {
    if (!is.null(state_dict)) {
        cells = tree$tip.label %>% intersect(names(state_dict))

        tree = tree %>% keep.tip(cells)
        state_dict = state_dict[cells]

        # Count categories in state_dict and filter by min_count
        state_counts <- table(state_dict)
        keep_states <- names(state_counts)[state_counts >= min_count]
        cells_keep <- names(state_dict)[state_dict %in% keep_states]
        tree <- keep.tip(tree, cells_keep)
        state_dict <- state_dict[cells_keep]

        # make sure X has the same cell ordering as W
        tree$edge.length = rep(1, nrow(tree$edge))
        tree$states = unname(state_dict[tree$tip.label])
        X <- PATH::catMat(tree$states)
    } else if (!is.null(feature_mat)) {
        cells = tree$tip.label %>% intersect(colnames(feature_mat))
        tree = tree %>% keep.tip(cells)
        # make sure X has the same cell ordering as W
        feature_mat = feature_mat[, tree$tip.label]
        X = t(feature_mat)
    } else {
        stop("Either state_dict or feature_mat must be provided.")
    }
    
    if (model == "sibling") {
        W = one_node_tree_dist2(tree, norm = norm)
    } else if (model == "bm") {
        W = bm_node_dist(tree, norm = norm)
    } else if (model == "inv") {
        W = PATH::inv_tree_dist(tree, norm = norm)
    } else {
        stop("Invalid model.")
    }

    tree$xcor <- xcor2(X, W)
    # tree$xcor <- PATH::xcor(X, W)

    res = tree$xcor$Z.score %>%
        reshape2::melt() %>%
        setNames(c('state_x', 'state_y', 'z')) %>%
        left_join(
            reshape2::melt(tree$xcor$phy_cor) %>%
                setNames(c('state_x', 'state_y', 'cor')),
            by = join_by(state_x, state_y)
        ) %>%
        left_join(
            reshape2::melt(tree$xcor$CI.lower) %>%
                setNames(c('state_x', 'state_y', 'lower')),
            by = join_by(state_x, state_y)
        ) %>%
        left_join(
            reshape2::melt(tree$xcor$CI.upper) %>%
                setNames(c('state_x', 'state_y', 'upper')),
            by = join_by(state_x, state_y)
        ) %>%
        mutate(p = 2*pnorm(-abs(z)))
    
    if (diag_only) {
        res = res %>% filter(state_x == state_y) %>% select(-state_y) %>% rename(state = state_x)
    }

    res = res %>% mutate(q = p.adjust(p, method = 'BH'))

    return(res)
    
}

run_path_inf = function(tree, state_dict, min_count = 0) {
    cells = tree$tip.label %>% intersect(names(state_dict))

    tree = tree %>% keep.tip(cells)
    state_dict = state_dict[cells]

    # Count categories in state_dict and filter by min_count
    state_counts <- table(state_dict)
    keep_states <- names(state_counts)[state_counts >= min_count]
    cells_keep <- names(state_dict)[state_dict %in% keep_states]
    tree <- keep.tip(tree, cells_keep)
    state_dict <- state_dict[cells_keep]

    # make sure X has the same cell ordering as W
    tree$edge.length = rep(1, nrow(tree$edge))
    tree$states = unname(state_dict[tree$tip.label])

    PATH:::PATH_inf(tree)
}

xcor2 <- function(data, weight.matrix) {
	#	data format: one row per cell and one column per measurement (e.g. one column per gene expression score)
	#	The ij'th element of the weight.matrix should reflect how closely related cells i and j are, 
	#	with larger numbers corresponding to closer relationships. The diagonal elements should be 0.

	t <- Matrix::t

	W <- methods::as(weight.matrix, "sparseMatrix")
	n <- nrow(data)
	
	d0 <- apply(as.matrix(data), 2, function(xx) xx - mean(xx))
	
	d1 <- (t(d0)%*%d0)/n
	d1.2 <- diag(d1)%*%t(diag(d1))
	d2 <- (t(d0^2)%*%(d0^2))/n
	
	x <- rep(1, n)
	w <- as.numeric(t(x)%*%W%*%x)
	S3 <- as.numeric(t(x)%*%(W * t(W))%*%x)
	S4 <- as.numeric(t(x)%*%(W * W)%*%x)
	S5 <- as.numeric(t(x)%*%W%*%W%*%x)
	S6 <- as.numeric(t(W%*%x)%*%(W%*%x) + t(t(W)%*%x)%*%(t(W)%*%x))
	S1 <- S4 + S3
	S2 <- 2*S5 + S6
	
	I <- as.matrix((t(d0)%*%W%*%d0)/(w*sqrt(d1.2)))
	
	E.I2 <- (-qlcMatrix::corSparse(data)/(n-1))
	
	V00 <- as.matrix((  ((d1^2 * n)/ (d1.2) ) * (2*(w^2 - S2 + S1) + (2*S3 - 2*S5)*(n-3) +
					       S3*(n-2)*(n-3) )  +  (( -d2 )/( d1.2 ))*(6*(w^2 - S2 + S1) +
				       (4*S1-2*S2)*(n-3) + S1*(n-2)*(n-3)  ) 
                      + n*((w^2 - S2 + S1) + (2*S4 -S6)*(n-3) + S4*(n-2)*(n-3) )  ) /
	((n-1)*(n-2)*(n-3)*(w^2) ) - 
	(  ( ( d1^2 )/( d1.2 ) ) * ( 1/((n-1)^2)  )  ))
	
	Z0 <- (I - E.I2)/sqrt(V00)
	pval <- stats::pnorm(-abs(Z0))*2

	#	Added: CI and SE (guard against tiny negative variance from numerical error)
	se <- sqrt(pmax(V00, 0))
	zcrit <- stats::qnorm(0.975)
	ci.lower <- I - zcrit * se
	ci.upper <- I + zcrit * se
	
	return(list(
		"phy_cor" = I,
		"Z.score" = Z0,
		"one.sided.pvalue" = pval,	# actually two-sided; kept name for backward compatibility
		"CI.lower" = ci.lower,
		"CI.upper" = ci.upper,
		"StdErr" = se
	))
}

##' Build a Brownian-motion style node distance matrix for a tree
##'
##' Compute a weighting matrix W suitable for PATH-style correlation tests.
##' This uses the phylogenetic variance-covariance matrix computed by
##' `ape::vcv` (after setting all edge lengths to 1), zeroes the diagonal,
##' optionally applies row-normalization (`PATH:::rowNorm`), and rescales
##' the matrix so entries sum to 1.
##'
##' @param tree A `phylo` object. Edge lengths will be set to 1 for this calculation.
##' @param norm Logical; if TRUE normalize rows with `PATH:::rowNorm` before scaling. Default: TRUE.
##' @return A numeric matrix (tips × tips) with row/column names matching `tree$tip.label`.
##' @examples
##' # W <- bm_node_dist(tree)
##' @export
bm_node_dist = function(tree, norm = TRUE) {
    tree$edge.length = rep(1, nrow(tree$edge))
    W = ape::vcv(tree)
    diag(W) = 0
    if (norm) {
        W = PATH:::rowNorm(W)
    }
    W <- W/sum(W)
    W
}

##' One-node sibling distance matrix
##'
##' Build a sparse tip-by-tip weighting matrix that connects tips sharing
##' the same immediate parent (i.e., sibling tips under the same internal
##' node). This wraps `one_node_depth2()` which returns a list including
##' a sparse matrix `W`; the returned matrix is optionally normalized and
##' scaled to sum to 1.
##'
##' @param tree A `phylo` object.
##' @param norm Logical; if TRUE apply `PATH:::rowNorm` before scaling. Default: TRUE.
##' @param exclude_root Logical; if TRUE internal nodes that are the root are excluded from sibling grouping. Default: FALSE.
##' @return A sparse numeric matrix (tips × tips) with row/column names matching `tree$tip.label`.
##' @examples
##' # W <- one_node_tree_dist2(tree)
##' @export
one_node_tree_dist2 = function(tree, norm = TRUE, exclude_root = FALSE) {
    w <- one_node_depth2(tree, exclude_root = exclude_root)$W
    if (norm == TRUE) {
        w <- PATH:::rowNorm(w)
    }
    w <- w/sum(w)
    rownames(w) = tree$tip.label
    colnames(w) = tree$tip.label
    return(w)
}

##' Extract one-node sibling adjacency and edge-pat statistics
##'
##' Scans a rooted phylogenetic tree for internal nodes that have two or
##' more tip descendants and builds a sparse tip×tip adjacency matrix
##' connecting sibling tip pairs. Also returns the mean parent-edge length
##' (mean.pat) and the vector of parent-edge lengths used (pats).
##'
##' @param tree A `phylo` object. If `tree$edge.length` is NULL it will be set to 1 for all edges.
##' @param exclude_root Logical; if TRUE, internal nodes equal to the root are excluded from grouping. Default: FALSE.
##' @return A list with elements:
##' * `W`: sparse symmetric adjacency matrix (tips × tips) connecting sibling tips
##' * `mean.pat`: mean parent-edge length times two
##' * `pats`: vector of parent-edge lengths for the included internal nodes
##' @examples
##' # out <- one_node_depth2(tree)
##' # W <- out$W; mean_pat <- out$mean.pat
##' @export
one_node_depth2 <- function(tree, exclude_root = FALSE) {

    if (is.null(tree$edge.length)) {
        tree$edge.length = rep(1, nrow(tree$edge))
    }

    bifur  <- castor::is_bifurcating(tree)
    Ntips  <- length(tree$tip.label)
    root   <- Ntips + 1L

    # find all edges that end in a tip
    tipnodes <- which(tree$edge[,2] %in% seq_len(Ntips))
    edges    <- tree$edge[tipnodes, , drop=FALSE]

    df <- data.frame(
        node   = edges[,1],
        tip    = edges[,2],
        pat    = tree$edge.length[tipnodes]
    )

    if (exclude_root) {
        df <- df %>% dplyr::filter(node != root)
    }

    # now group by internal node and keep only those with ≥2 tips
    df2 <- df %>%
        dplyr::group_by(node) %>%
        dplyr::mutate(clades = list(tip)) %>%
        dplyr::rowwise() %>%
        dplyr::filter(length(clades) >= 2)

    # build the tip-tip adjacency
    indices <- df2 %>%
        dplyr::select(-tip, -pat) %>%
        dplyr::distinct() %>%
        dplyr::pull(clades) %>%
    {
        if (!bifur) {
            lapply(., function(x) t(utils::combn(x, 2))) %>%
                do.call(rbind, .)
        } else {
            do.call(rbind, .)
        }
    }

    mat <- Matrix::sparseMatrix(
        i = c(indices[,1], indices[,2]),
        j = c(indices[,2], indices[,1]),
        dims = c(Ntips, Ntips)
    )

    mp <- mean(df2$pat) * 2

    list(W = mat, mean.pat = mp, pats = df2$pat)
}


#' Run phylogenetic signal analysis on a trait matrix
#'
#' This function computes the Cmean statistic (and its p-value) for each feature
#' in a trait matrix against a given phylogenetic tree. It handles intersection
#' of cells, removal of zero-variance features, and parallel execution.
#'
#' @param trait_mat A numeric matrix of traits/features (rows) x cells (columns).
#' @param tree A `phylo` object representing the cell lineage tree.
#' @param reps Integer; number of permutations for the randomization test. Default: 999.
#' @param ncores Integer; number of cores to use for parallel processing. Default: 1.
#' @param verbose Logical; whether to print progress messages. Default: TRUE.
#' @return A data frame with columns:
#'   - `gene`: Feature name
#'   - `cmean`: The Cmean statistic
#'   - `pval`: P-value from the randomization test
#'   - `qval`: BH-adjusted p-value
#' @export
run_signal <- function(trait_mat, tree, reps = 999, ncores = 1, verbose = TRUE) {

    RhpcBLASctl::blas_set_num_threads(1)

    trait_mat <- as.matrix(trait_mat)

    cells_common = intersect(colnames(trait_mat), tree$tip.label)
    tree = tree %>% ape::keep.tip(cells_common)
    tree$node.label = NULL

    trait_mat = trait_mat[,cells_common,drop=FALSE]

    feature_variances <- apply(trait_mat, 1, var, na.rm = TRUE)
    zero_var_features <- which(feature_variances == 0)

    if (length(zero_var_features) > 0) {
        message(sprintf("Removing %d features with zero variance", length(zero_var_features)))
        trait_mat <- trait_mat[-zero_var_features, , drop = FALSE]
    }

    trait_dat = t(trait_mat)
    
    df_test <- mclapply(
        colnames(trait_dat),
        mc.cores = ncores,
        function(trait) {
            if (verbose) {
                message(trait)
            }

            p4d_i <- phylobase::phylo4d(tree, trait_dat[, trait, drop = FALSE]) %>% suppressMessages()
            test_dat = phylosignal::phyloSignal(p4d = p4d_i,
                    method = "Cmean",
                    reps   = reps)
            df_test = full_join(
                    test_dat$stat %>% tibble::rownames_to_column('gene') %>% rename(cmean = Cmean),
                    test_dat$pvalue %>% tibble::rownames_to_column('gene') %>% rename(pval = Cmean),
                    by = join_by(gene)
                )
        }) %>%
        bind_rows() %>%
        mutate(gene = colnames(trait_dat)) %>%
        mutate(qval = p.adjust(pval, method = 'BH')) %>%
        arrange(pval)
    
    return(df_test)
}


#' Plot lineage coupling heatmap
#'
#' Visualizes the pairwise correlation or Z-scores between states as a heatmap.
#' Optionally performs hierarchical clustering to order states and can display
#' the clustering dendrogram.
#'
#' @param res_path Data frame containing results (e.g. from `run_path`). Must contain columns `state_x`, `state_y`, `cor`, `z`, `q`.
#' @param state_order Optional character vector defining the order of states. If NULL, hierarchical clustering is performed.
#' @param limits Numeric vector of length 2 defining the color scale limits. Default: c(-0.5, 0.5).
#' @param statistic Character; choose whether to display "cor" or "z" from `res_path`. Default: "cor".
#' @param mark_signif Logical; if TRUE, overlay points for high Z-scores (|z|>2.5) and tiles for significant q-values. Default: TRUE.
#' @param signif_only Logical; if TRUE, mask non-significant cells (q > q_thres). Default: FALSE.
#' @param q_thres Numeric; threshold for q-value significance. Default: 0.05.
#' @param rev_order Logical; if TRUE and state_order is NULL, reverse the clustered order. Default: FALSE.
#' @param highlight Optional character vector of state names to highlight in axis labels.
#' @param highlight_color Color string for highlighted labels. Default: "firebrick".
#' @param plot_tree Logical; if TRUE, plot the dendrogram above the heatmap (requires `ggtree` and `patchwork`). Default: FALSE.
#' @param dend_height Numeric; relative height of the dendrogram. Default: 0.3.
#' @return A ggplot object (or patchwork object if `plot_tree` is TRUE).
#' @export
plot_res_path = function(
    res_path, state_order = NULL, limits = c(-0.5, 0.5), statistic = c("cor", "z"), mark_signif = TRUE, 
    signif_only = FALSE, q_thres = 0.05, rev_order = FALSE, highlight = NULL, highlight_color = "firebrick",
     plot_tree = FALSE, dend_height = 0.3
) {

    statistic <- match.arg(statistic)

    res_path = res_path %>% filter(!is.na(cor))

    if (is.null(state_order)) {
        cor_mat = res_path %>%
            reshape2::dcast(state_x ~ state_y, value.var = 'cor') %>% 
            tibble::column_to_rownames('state_x')
        hc <- hclust(as.dist(1-cor_mat), method = "complete")
        state_order = hc$labels[hc$order]
        if (rev_order) {
            state_order = rev(state_order)
        }
    }

    if (statistic == "z") {
        if (!"z" %in% colnames(res_path)) {
            stop("`res_path` must contain a `z` column to plot statistic = 'z'.")
        }
        res_path <- res_path %>% mutate(plot_val = z)
        guide_title <- "Coupling\nZ-score"
    } else {
        res_path <- res_path %>% mutate(plot_val = cor)
        guide_title <- "Coupling\ncorrelation"
    }
    
    p = res_path %>%
        filter(state_x %in% state_order) %>%
        filter(state_y %in% state_order) %>%
        mutate(state_x = factor(state_x, state_order)) %>%
        mutate(state_y = factor(state_y, state_order)) %>%
        mutate(cor_plot = ifelse(q > q_thres & signif_only, NA_real_, plot_val)) %>%
        ggplot(
            aes(x = state_x, y = state_y, fill = cor_plot)
        ) +
        geom_raster() +
        scale_fill_gradient2(
            low = 'blue', high = 'red', midpoint = 0,
            limits = limits, oob = scales::oob_squish, na.value = 'white'
        ) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        scale_y_discrete(position = "right") +
        guides(fill = guide_colorbar(title = guide_title))

    # If highlighted features are provided, color their tick labels instead of dropping others
    if (!is.null(highlight)) {
        hl <- intersect(highlight, state_order)
        if (length(hl) > 0) {
            if (requireNamespace("ggtext", quietly = TRUE)) {
                lab_fun <- function(x) ifelse(x %in% hl, paste0("<span style='color:", highlight_color, "; font-weight:700;'>", x, "</span>"), x)
                p <- p +
                    scale_x_discrete(labels = lab_fun, drop = FALSE) +
                    scale_y_discrete(labels = lab_fun, position = "right", drop = FALSE) +
                    theme(axis.text.x = ggtext::element_markdown(angle = 30, hjust = 1),
                          axis.text.y = ggtext::element_markdown())
            } else {
                # fallback: bold highlighted labels using plotmath
                lab_expr <- function(x) {
                    sapply(x, function(s) if (s %in% hl) paste0('bold("', s, '")') else paste0('"', s, '"'))
                }
                p <- p +
                    scale_x_discrete(labels = function(x) parse(text = lab_expr(x)), drop = FALSE) +
                    scale_y_discrete(labels = function(y) parse(text = lab_expr(y)), position = "right", drop = FALSE)
            }
        }
    }
    
    if (mark_signif) {
        p = p +
        geom_point(
            data = . %>% filter(abs(z) > 2.5),
            inherit.aes = F,
            aes(x = state_x, y = state_y)
        ) +
        geom_tile(
            data = . %>% filter(q < q_thres),
            inherit.aes = F,
            color = 'black',
            fill = NA,
            linewidth = 0.25,
            aes(x = state_x, y = state_y)
        )
    }

    # Optionally plot the hierarchical clustering tree above the heatmap
    if (isTRUE(plot_tree)) {

        # convert hclust to phylo and plot with ggtree
        tree_phy <- ape::as.phylo(hc)
        dend_plot <- ggtree::ggtree(tree_phy, branch.length = 'none') +
            theme_tree2() +
            theme(
                plot.margin = margin(0,0,0,0),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()
            )

        # combine dendrogram above heatmap using patchwork
        combined <- patchwork::wrap_plots(dend_plot, p, nrow = 1, widths = c(dend_height, 1))
        return(combined)
    }

    return(p)
}

#' Plot lineage coupling triangular heatmap
#'
#' Visualizes the pairwise correlation or Z-scores between states as a triangular
#' (diamond-shaped) heatmap. This is useful for symmetric matrices.
#'
#' @param res_path Data frame containing results. Must contain columns `state_x`, `state_y`, `cor`, `z`, `q`.
#' @param state_order Optional character vector defining the order of states. If NULL, hierarchical clustering is performed.
#' @param limits Numeric vector of length 2 defining the color scale limits. Default: c(-0.5, 0.5).
#' @param statistic Character; choose whether to display "cor" or "z" from `res_path`. Default: "cor".
#' @param mark_signif Logical; if TRUE, overlay points for high Z-scores and outlines for significant q-values. Default: TRUE.
#' @param signif_only Logical; if TRUE, mask non-significant cells (q > q_thres). Default: FALSE.
#' @param q_thres Numeric; threshold for q-value significance. Default: 0.05.
#' @param rev_order Logical; if TRUE and state_order is NULL, reverse the clustered order. Default: FALSE.
#' @param highlight Optional character vector of state names to highlight in axis labels.
#' @param highlight_color Color string for highlighted labels. Default: "firebrick".
#' @return A ggplot object.
#' @export

plot_coupling_triangle = function(
    res_path, state_order = NULL, limits = c(-0.5, 0.5), statistic = c("cor", "z"), mark_signif = TRUE,
    signif_only = FALSE, q_thres = 0.05, rev_order = FALSE, highlight = NULL, highlight_color = "firebrick"
) {

    statistic <- match.arg(statistic)
    res_path = res_path %>% filter(!is.na(cor))

    if (is.null(state_order)) {
        cor_mat = res_path %>%
            reshape2::dcast(state_x ~ state_y, value.var = 'cor') %>% 
            tibble::column_to_rownames('state_x')
        cor_mat[is.na(cor_mat)] = 0
        hc <- hclust(as.dist(1-cor_mat), method = "complete")
        state_order = hc$labels[hc$order]
        if (rev_order) {
            state_order = rev(state_order)
        }
    }

    if (statistic == "z") {
        if (!"z" %in% colnames(res_path)) {
            stop("`res_path` must contain a `z` column to plot statistic = 'z'.")
        }
        res_path <- res_path %>% mutate(plot_val = z)
        guide_title <- "Coupling\nZ-score"
    } else {
        res_path <- res_path %>% mutate(plot_val = cor)
        guide_title <- "Coupling\ncorrelation"
    }
    
    # Map states to indices
    state_map = setNames(seq_along(state_order), state_order)
    
    df_plot = res_path %>%
        filter(state_x %in% state_order, state_y %in% state_order) %>%
        mutate(
            idx_x = state_map[as.character(state_x)],
            idx_y = state_map[as.character(state_y)]
        ) %>%
        # Ensure idx_x <= idx_y (upper triangle)
        rowwise() %>%
        mutate(
            i = min(idx_x, idx_y),
            j = max(idx_x, idx_y)
        ) %>%
        ungroup() %>%
        # Remove duplicates if input had both (A,B) and (B,A)
        distinct(i, j, .keep_all = TRUE) %>%
        mutate(cor_plot = ifelse(q > q_thres & signif_only, NA_real_, plot_val))

    # Expand to polygon coordinates (diamonds)
    df_poly = df_plot[rep(1:nrow(df_plot), each = 4), ]
    df_poly$vertex = rep(1:4, nrow(df_plot))
    
    df_poly = df_poly %>%
        mutate(
            x = case_when(
                vertex == 1 ~ (i + j - 1)/2,
                vertex == 2 ~ (i + j)/2,
                vertex == 3 ~ (i + j + 1)/2,
                vertex == 4 ~ (i + j)/2
            ),
            y = case_when(
                vertex == 1 ~ (j - i)/2,
                vertex == 2 ~ (j - i - 1)/2,
                vertex == 3 ~ (j - i)/2,
                vertex == 4 ~ (j - i + 1)/2
            ),
            group_id = paste(i, j, sep = "-")
        )

    # identify boundary edges (appear only once across all diamonds)
    boundary_edges <- df_poly %>%
        arrange(group_id, vertex) %>%
        group_by(group_id) %>%
        mutate(
            x_next = lead(x, default = first(x)),
            y_next = lead(y, default = first(y))
        ) %>%
        ungroup() %>%
        mutate(
            edge_key = paste(
                sprintf("%.6f", pmin(x, x_next)),
                sprintf("%.6f", pmax(x, x_next)),
                sprintf("%.6f", pmin(y, y_next)),
                sprintf("%.6f", pmax(y, y_next)),
                sep = "|"
            )
        ) %>%
        add_count(edge_key, name = "edge_n") %>%
        filter(edge_n == 1)

    p = ggplot(df_poly, aes(x = x, y = y, group = group_id, fill = cor_plot)) +
        geom_polygon(color = "black", linewidth = 0.15) +
        scale_fill_gradient2(
            low = 'blue', high = 'red', mid = 'white', midpoint = 0,
            limits = limits, oob = scales::oob_squish, na.value = 'white'
        ) +
        geom_segment(
            data = boundary_edges,
            aes(x = x, y = y, xend = x_next, yend = y_next),
            inherit.aes = FALSE,
            linewidth = 0.4,
            color = "black",
            lineend = "round"
        ) +
        theme_void() +
        theme(
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            legend.position = "right"
        ) +
        guides(fill = guide_colorbar(title = guide_title)) +
        scale_x_continuous(
            breaks = seq_along(state_order),
            labels = state_order
        ) +
        coord_fixed()

    # Highlight logic
    if (!is.null(highlight)) {
        hl <- intersect(highlight, state_order)
        if (length(hl) > 0) {
            if (requireNamespace("ggtext", quietly = TRUE)) {
                lab_fun <- function(x) ifelse(x %in% hl, paste0("<span style='color:", highlight_color, "; font-weight:700;'>", x, "</span>"), x)
                p <- p + theme(axis.text.x = ggtext::element_markdown())
                p <- p + scale_x_continuous(breaks = seq_along(state_order), labels = lab_fun(state_order))
            } else {
                # fallback: bold highlighted labels using plotmath
                lab_expr <- function(x) {
                    sapply(x, function(s) if (s %in% hl) paste0('bold("', s, '")') else paste0('"', s, '"'))
                }
                p <- p + scale_x_continuous(breaks = seq_along(state_order), labels = function(x) parse(text = lab_expr(x)))
            }
        }
    }
    
    if (mark_signif) {
        # Centers for points
        df_centers = df_plot %>%
            mutate(
                x = (i + j) / 2,
                y = (j - i) / 2
            )
            
        p = p +
            geom_point(
                data = df_centers %>% filter(q < q_thres),
                aes(x = x, y = y, group = NULL, fill = NULL),
                inherit.aes = FALSE,
                size = 1,
                pch = 8,
                stroke = 0.25
            )
    }

    return(p)
}

#' Plot self-coupling estimates with confidence intervals
#'
#' Extracts the diagonal entries from a PATH result table and visualizes each
#' state's self-coupling estimate (`cor`) together with its lower/upper
#' confidence interval using `geom_pointrange`.
#'
#' @param res_path Data frame produced by `run_path` (or equivalent) that must
#'   contain columns `state_x`, `state_y`, `cor`, `lower`, and `upper`.
#' @param state_order Optional character vector specifying the plotting order of
#'   states. When `NULL`, states are sorted alphabetically.
#' @param q_thres Numeric significance cutoff; states with diagonal q-values
#'   below this threshold are highlighted automatically. Default: 0.05.
#' @param highlight_color Color applied to highlighted states. Default:
#'   "firebrick".
#' @return A ggplot object showing the point estimates with confidence ranges.
#' @export
plot_self_coupling <- function(
    res_path,
    state_order = NULL,
    q_thres = 0.05,
    dot_size = 0.5,
    line_width = 0.5,
    highlight_color = "firebrick"
) {

    required_cols <- c("state_x", "state_y", "cor", "lower", "upper", "q")
    missing_cols <- setdiff(required_cols, colnames(res_path))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "res_path is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    diag_df <- res_path %>%
        filter(state_x == state_y) %>%
        transmute(
            state = state_x,
            cor = cor,
            lower = lower,
            upper = upper,
            q = q
        ) %>%
        distinct(state, .keep_all = TRUE)

    if (nrow(diag_df) == 0) {
        stop("No diagonal entries (state_x == state_y) were found in res_path.")
    }

    if (is.null(state_order)) {
        state_order <- sort(unique(diag_df$state))
    } else {
        # keep any states missing from the supplied order appended to the end
        missing_states <- setdiff(diag_df$state, state_order)
        state_order <- c(state_order, missing_states)
    }

    diag_df <- diag_df %>%
        mutate(state = factor(state, levels = state_order)) %>%
        arrange(state) %>%
        mutate(idx = as.numeric(state))

    n_states <- length(state_order)
    y_min_raw <- min(diag_df$cor, na.rm = TRUE)
    y_max_raw <- max(diag_df$cor, na.rm = TRUE)

    y_range <- y_max_raw - y_min_raw
    pad <- 0.2 * y_range

    y_min <- 0
    y_max <- y_max_raw + pad
    border_df <- tibble::tibble(
        x = c(0.5, n_states + 0.5, n_states + 0.5, 0.5, 0.5),
        y = c(y_min, y_min, y_max, y_max, y_min)
    )

    baseline_df <- tibble::tibble(
        x = 0.5,
        xend = n_states + 0.5,
        y = 0,
        yend = 0
    )

    diag_df <- diag_df %>%
        mutate(color_group = ifelse(!is.na(q) & q < q_thres, "highlight", "default"))

    color_values <- c(default = "grey70")
    if (any(diag_df$color_group == "highlight")) {
        color_values <- c(color_values, highlight = highlight_color)
    }

    p <- ggplot(diag_df, aes(x = idx, y = cor)) +
        geom_segment(
            data = baseline_df,
            aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE,
            color = "grey80",
            linewidth = 0.25
        ) +
        geom_col(
            aes(fill = color_group),
            width = 0.8,
            color = "black",
            linewidth = 0.2
        ) +
        geom_text(
            data = diag_df %>% filter(color_group == "highlight"),
            aes(x = idx, y = cor, label = "*"),
            inherit.aes = FALSE,
            vjust = 0.25,
            size = 3,
            fontface = "bold"
        ) +
        scale_x_continuous(
            breaks = seq_along(state_order),
            labels = state_order,
            expand = expansion(add = c(0, 0))
        ) +
        # scale_y_continuous(
        #     limits = c(y_min, y_max),
        #     sec.axis = dup_axis(
        #         name = "Heritability",
        #         labels = NULL,
        #         breaks = NULL
        #     ),
        #     expand = expansion(add = c(0, 0))
        # ) +
        scale_fill_manual(values = color_values, guide = "none") +
        labs(
            x = NULL,
            y = NULL
        ) +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 2)),
            axis.text.y = element_text(size = 6),
            plot.margin = margin(t = 5, r = 10, b = 10, l = 50),
            axis.title.y = element_blank(),
            axis.title.y.right = element_text(angle = 270, vjust = 0.5, hjust = 0.5, size = 9, margin = margin(l = 4)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_blank()
        ) +
        geom_path(
            data = border_df,
            aes(x = x, y = y),
            inherit.aes = FALSE,
            linewidth = 0.3,
            color = "black"
        )

    return(p)
}


plot_coupling_ribbons <- function(
	A,
    state_order = NULL,
	alpha = 0.7,
	edge_color = "grey30",   # kept but not used anymore
	node_color = "black",
	linewidth_range = c(0.5, 6),
	node_size_range = c(2, 6),
    min_prob = 0,
    node_fill_map = NULL,
    node_size_map = NULL
) {
    if (is.null(state_order)) {
        state_order <- colnames(A)
    }

	df_trans = as.data.frame(A) %>%
		tibble::rownames_to_column("state_x") %>%
		reshape2::melt(id.var = "state_x", variable.name = "state_y", value.name = "prob")

	df <- df_trans %>%
		mutate(
			state_x = factor(as.character(state_x), levels = state_order),
			state_y = factor(as.character(state_y), levels = state_order)
		) %>%
		filter(!is.na(state_x), !is.na(state_y)) %>%
		mutate(
			idx_x = as.numeric(state_x),
			idx_y = as.numeric(state_y)
		)

	# self-coupling for node size (unchanged)
	self_df <- df %>%
		filter(idx_x == idx_y) %>%
		group_by(idx_x) %>%
		summarise(self_prob = sum(prob, na.rm = TRUE), .groups = "drop")

	# combine unordered pairs, but keep single prob (assume symmetric)
	pair_df <- df %>%
		filter(idx_x != idx_y) %>%
		mutate(
			idx_left = pmin(idx_x, idx_y),
			idx_right = pmax(idx_x, idx_y),
			pair = paste(idx_left, idx_right, sep = "||")
		) %>%
		group_by(pair, idx_left, idx_right) %>%
		summarise(weight = sum(prob), .groups = "drop") %>%
		filter(weight > min_prob) %>%
		mutate(
			x = idx_left,
			xend = idx_right
		)

	# --- NEW: stretch weights + map to linewidth, color, and alpha ----
	if (nrow(pair_df) > 0) {
		# optional: non-linear emphasis on strong edges
		w_raw <- pair_df$weight
		w_emph <- w_raw^1.5  # gamma > 1 emphasizes big edges; adjust if you want

		pair_df <- pair_df %>%
			mutate(
				w_plot  = scales::rescale(w_emph, to = linewidth_range),
				col_val = scales::rescale(w_raw, to = c(0, 1)),
				alpha_val = scales::rescale(w_raw, to = c(0.15, 0.9))
            ) %>%
            arrange(weight)  # draw weakest first so strong ribbons stay on top
	}

    state_df <- tibble::tibble(
            x = seq_along(state_order),
            y = 0,
            state = state_order
        ) %>%
		left_join(self_df, by = c("x" = "idx_x"))

    # optional overrides for node size
    state_df$size <- mean(node_size_range)
    if (!is.null(node_size_map)) {
        state_df$size <- node_size_map[state_df$state]
    }

    # determine fill colors per node
    state_df$fill <- "white"
    if (!is.null(node_fill_map)) {
        state_df$fill <- node_fill_map[state_df$state]
    }

	g <- ggplot() +
		geom_curve(
			data = pair_df,
			aes(
				x = x,
				y = 0,
				xend = xend,
				yend = 0,
				linewidth = w_plot,
				color    = col_val,
				alpha    = alpha_val
			),
			curvature = 0.3,
			lineend = "round"
		) +
        geom_point(
            data = state_df,
            aes(x = x, y = y, size = size, fill = fill),
            shape = 21,
            stroke = 0.4,
            color = node_color
        ) +
		scale_size_identity() +
        scale_fill_identity(guide = "none") +
		scale_linewidth_identity(
			guide = "none"
		) +
		scale_color_gradient(
			name = "Transition\nstrength",
			low  = "#deebf7",   # light blue (ColorBrewer Blues)
			high = "#08519c"   # deep blue
		) +
		scale_alpha_identity(guide = "none") +
		scale_x_continuous(
			breaks = seq_along(state_order),
			labels = state_order,
			expand = expansion(add = 0.5)
		) +
		# ylim(-0.6, 0.5) +
		theme_minimal(base_size = 11) +
		theme(
			axis.title = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks = element_blank(),
			panel.grid = element_blank(),
			axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
			legend.position = "right",
			legend.justification = "center"
		)

	return(g)
}

plot_marker_heatmap <- function(
    df_z,
    state_order,
    highlight = c(),
    x_limits = NULL,
    z_limits = c(-1,1),
    x_breaks = NULL,
    label_size = 2,
    show_xlabel = TRUE,
    y_title = NULL,
    left_pad = 0.35,
    legend_title = "Normalized\nexpression",
    plot_title = NULL
) {

    # make sure cluster factor order matches p_coupling
    df_plot <- df_z %>%
        mutate(
            cluster = as.integer(factor(cluster, levels = state_order))
        )

    pad_val <- if (length(highlight) > 0) left_pad else 0

    if (is.null(x_limits)) {
        x_limits <- c(0.5 - pad_val, length(state_order) + 0.5)
    }

    if (is.null(x_breaks)) {
        x_breaks <- seq_along(state_order)
    }

    label_anchor <- x_limits[1] + pad_val

    # one row per highlighted gene, anchored just inside the padded left margin
    label_df <- df_plot %>%
        filter(feature %in% highlight) %>%
        distinct(feature) %>%
        mutate(
            cluster = label_anchor
        )

    axis_title_right <- if (!is.null(y_title) && nzchar(y_title)) {
        element_text(angle = 270, vjust = 0.5, hjust = 0.5, size = 9, margin = margin(l = 4))
    } else {
        element_blank()
    }

    n_states <- length(state_order)
    feature_levels <- levels(df_plot$feature)
    if (is.null(feature_levels)) {
        feature_levels <- unique(df_plot$feature)
    }
    n_features <- length(feature_levels)

    p_markers_rna <- ggplot(df_plot, aes(x = cluster, y = feature, fill = score)) +
        geom_tile() +
        scale_fill_gradient2(
            low = "#018571", high = "#A01F7A", midpoint = 0,
            limits = z_limits, oob = scales::oob_squish
        ) +
        scale_x_continuous(breaks = x_breaks, labels = state_order, expand = c(0, 0)) +
        scale_y_discrete(position = "right", expand = c(0, 0)) +
        coord_cartesian(xlim = x_limits, clip = "off") +
        theme_void() +
        theme(
            axis.text.x = if (show_xlabel) element_text(angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 2)) else element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.y.right = axis_title_right,
            plot.margin = margin(t = 0, r = 10, b = 10, l = 50),
            plot.title = element_text(hjust = 0.5, size = 12)
        ) +
        labs(y = y_title, title = plot_title) +
        geom_text_repel(
            data = label_df,
            aes(x = cluster, y = feature, label = feature),
            hjust = 1,
            nudge_x = -0.2,
            xlim = c(-Inf, x_limits[1]),
            direction = "y",
            segment.inflect = TRUE,
            segment.color = "black",
            segment.size = 0.1,
            segment.curvature = 0.2,
            segment.angle = 90,
            segment.ncp = 1,
            segment.square = TRUE,
            min.segment.length = 0,
            box.padding = 0,
            point.padding = 0,
            size = label_size,
            inherit.aes = FALSE
        ) +
        guides(fill = guide_colorbar(title = legend_title))

    if (n_features > 0) {
        border_df <- tibble::tibble(
            x = c(0.5, n_states + 0.5, n_states + 0.5, 0.5, 0.5),
            y = c(0.5, 0.5, n_features + 0.5, n_features + 0.5, 0.5)
        )

        p_markers_rna <- p_markers_rna +
            geom_path(
                data = border_df,
                aes(x = x, y = y),
                inherit.aes = FALSE,
                linewidth = 0.3,
                color = "black"
            )
    }

    return(p_markers_rna)
}


#' Plot PATH coupling panel
#'
#' Assemble the triangular PATH heatmap, ribbon diagram, and optional marker
#' heatmaps into a single stacked panel. Marker inputs now accept either a
#' single data frame (for one modality) or a list of data frames (to plot
#' multiple modalities in separate heatmaps).
#'
#' @param marker_df Either one data frame (for a single heatmap) or a list of
#'   data frames, each already formatted for `plot_marker_heatmap`.
#' @param marker_highlight Character vector of genes to label, or a list of
#'   character vectors matching the entries in `marker_df` when multiple panels
#'   are provided.
#' @param marker_height Numeric scalar or vector giving the relative height for
#'   each marker panel (recycled when a single value is supplied).
#' @param marker_show_xlabel Logical; if FALSE, suppress x-axis tick labels on
#'   the marker heatmaps (useful when stacking many panels). Default: TRUE.
#' @param marker_left_pad Numeric amount of left padding (in x-axis units)
#'   reserved when marker highlights are shown. Default: 0.35.
#' @param panel_spacing Numeric height (relative to panel heights) of blank space
#'   inserted between stacked panels. Default: 0 (no extra spacing).
plot_coupling_panel <- function(
    path_res,
    A,
    state_order = NULL,
    linewidth_range = c(0.5, 6),
    limits = c(-0.5, 0.5),
    self_coupling_limits = NULL,
    statistic = c("cor", "z"),
    node_fill_map = NULL,
    node_size_map = NULL,
    mark_signif = TRUE,
    marker_df = NULL,
    marker_highlight = character(),
    marker_height = 1.5,
    marker_label_size = 1,
    marker_show_xlabel = TRUE,
    marker_z_limits = c(-1, 1),
    marker_left_pad = 0.35,
    panel_spacing = 0,
    self_coupling_height = 0.9
) {

    if (is.null(state_order)) {
        state_order <- colnames(A)
    }

    statistic <- match.arg(statistic)

    has_highlights <- FALSE
    if (is.list(marker_highlight) && !is.data.frame(marker_highlight)) {
        has_highlights <- any(lengths(marker_highlight) > 0)
    } else if (!is.null(marker_highlight)) {
        has_highlights <- length(marker_highlight) > 0
    }

    left_pad <- if (has_highlights) marker_left_pad else 0
    xlim_common <- c(0.5 - left_pad, length(state_order) + 0.5)
    x_breaks <- seq_along(state_order)

    p1 <- plot_coupling_triangle(
            path_res,
            state_order = state_order,
            statistic = statistic,
            limits = limits,
            mark_signif = mark_signif
        ) +
        scale_x_continuous(
            breaks = x_breaks,
            labels = state_order,
            expand = c(0, 0)
        ) +
        coord_cartesian(xlim = xlim_common, clip = "off") +
        theme_void() +
        theme(
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
        )

    p2 <- plot_coupling_ribbons(
            A,
            state_order = state_order,
            linewidth_range = linewidth_range,
            node_fill_map = node_fill_map,
            node_size_map = node_size_map
        ) +
        scale_x_continuous(
            breaks = x_breaks,
            labels = state_order,
            expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(0, 0), expand = expansion(mult = c(0, 0))) +
        coord_cartesian(xlim = xlim_common, clip = "off") +
        theme_void() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(t = 0, r = 10, b = 30, l = 10)
        )

    panel_plots <- list()
    panel_heights <- numeric(0)

    add_panel <- function(panel, height) {
        panel_plots[[length(panel_plots) + 1]] <<- panel
        panel_heights <<- c(panel_heights, height)
    }

    add_panel(p1, 3)
    add_panel(p2, 0.32)

    p_self <- plot_self_coupling(
            path_res,
            state_order = state_order
        ) +
        coord_cartesian(xlim = xlim_common, ylim = self_coupling_limits, clip = "off") +
        theme(
            plot.margin = margin(t = 0, r = 10, b = 10, l = 10)
        )

    add_panel(p_self, self_coupling_height)

    if (!is.null(marker_df)) {
        if (is.data.frame(marker_df)) {
            marker_list <- list(marker_df)
        } else if (is.list(marker_df)) {
            marker_list <- marker_df
        } else {
            stop("marker_df must be a data frame or a list of data frames.")
        }

        marker_titles <- names(marker_list)

        highlight_list <- marker_highlight
        if (is.null(highlight_list)) {
            highlight_list <- list(character())
        } else if (is.list(highlight_list) && !is.data.frame(highlight_list)) {
            highlight_list <- highlight_list
        } else {
            highlight_list <- list(as.character(highlight_list))
        }

        if (length(highlight_list) == 1 && length(marker_list) > 1) {
            highlight_list <- rep(highlight_list, length(marker_list))
        } else if (length(highlight_list) != length(marker_list)) {
            stop("marker_highlight must be a list the same length as marker_df when providing multiple panels.")
        }

        marker_height_vec <- marker_height
        if (length(marker_height_vec) == 1) {
            marker_height_vec <- rep(marker_height_vec, length(marker_list))
        } else if (length(marker_height_vec) != length(marker_list)) {
            stop("marker_height must be length 1 or match the number of marker_df entries.")
        }

        marker_xlabel_vec <- marker_show_xlabel
        if (length(marker_xlabel_vec) == 1) {
            marker_xlabel_vec <- rep(marker_xlabel_vec, length(marker_list))
        } else if (length(marker_xlabel_vec) != length(marker_list)) {
            stop("marker_show_xlabel must be length 1 or match the number of marker_df entries.")
        }
        marker_xlabel_vec <- as.logical(marker_xlabel_vec)

        for (i in seq_along(marker_list)) {
            p_markers <- plot_marker_heatmap(
                df_z = marker_list[[i]],
                state_order = state_order,
                highlight = highlight_list[[i]],
                x_limits = xlim_common,
                x_breaks = x_breaks,
                label_size = marker_label_size,
                show_xlabel = marker_xlabel_vec[[i]],
                z_limits = marker_z_limits,
                y_title = if (!is.null(marker_titles) && nzchar(marker_titles[i])) marker_titles[i] else NULL,
                left_pad = marker_left_pad
            ) +
            coord_cartesian(xlim = xlim_common, clip = "off") +
            theme(
                plot.margin = margin(t = 0, r = 10, b = 10, l = 50)
            )

            add_panel(p_markers, marker_height_vec[[i]])
        }
    }

    plots <- panel_plots
    heights <- panel_heights

    if (panel_spacing > 0 && length(plots) > 1) {
        n_panels <- length(plots)
        new_plots <- vector("list", 2 * n_panels - 1)
        new_heights <- numeric(2 * n_panels - 1)
        for (i in seq_len(n_panels)) {
            idx <- 2 * i - 1
            new_plots[[idx]] <- plots[[i]]
            new_heights[idx] <- heights[i]
            if (i < n_panels) {
                new_plots[[idx + 1]] <- patchwork::plot_spacer()
                new_heights[idx + 1] <- panel_spacing
            }
        }
        plots <- new_plots
        heights <- new_heights
    }

    combined <- patchwork::wrap_plots(
        plots = plots,
        ncol = 1,
        heights = heights,
        guides = "collect"
    )

    combined
}

##' 1-step bidirectional proximity to target states
##'
##' Compute, for each state, its 1-step bidirectional connectivity to a set of
##' target states given a square transition matrix. Connectivity for state *i*
##' is based on the symmetrised weights (A[i, j] + A[j, i]) towards target
##' states. Optionally, scores are normalised by the total bidirectional weight
##' of each state.
##'
##' Target states are listed first in the output and ordered by their proximity
##' to non-target states. Non-target states are ordered by their proximity to
##' target states.
##'
##' @param A Square numeric matrix where `A[i, j]` is the transition
##'   probability or weight from state *i* to state *j*.
##' @param target Character, logical, or numeric vector indicating target
##'   states. If character, values must match `colnames(A)`. Logical/numeric
##'   vectors must have length `nrow(A)`; non-zero / TRUE entries are treated
##'   as targets.
##' @param normalize Logical; if TRUE (default), divide raw connectivity scores
##'   by the total bidirectional weight of each state. If FALSE, return
##'   unnormalised scores.
##'
##' @return A data frame with one row per state and columns:
##'   \describe{
##'     \item{state}{State name (from `rownames(A)` if available).}
##'     \item{score}{1-step proximity to target states (normalised if
##'       `normalize = TRUE`).}
##'     \item{score_nt}{1-step proximity to non-target states (normalised if
##'       `normalize = TRUE`).}
##'     \item{raw_score}{Unnormalised proximity to target states.}
##'     \item{raw_score_nt}{Unnormalised proximity to non-target states.}
##'     \item{is_target}{Logical flag indicating whether the state is a target.}
##'   }
##'   Target states appear first (ordered by increasing proximity to
##'   non-targets), followed by non-target states (ordered by decreasing
##'   proximity to targets).
##'
##' @export
target_proximity <- function(A, target, normalize = TRUE) {
	# A: square matrix, A[i, j] = transition prob/weight i -> j
	# target: character, logical, or numeric vector indicating target states

	A <- as.matrix(A)
	n <- nrow(A)
	if (ncol(A) != n) stop("A must be square")

	# build target indicator vector r (length n)
	if (is.character(target)) {
		if (is.null(colnames(A))) stop("A must have row/colnames when target is character")
		r <- as.numeric(colnames(A) %in% target)
	} else if (is.logical(target) || is.numeric(target)) {
		if (length(target) != n) stop("target vector must have length nrow(A)")
		r <- as.numeric(target != 0)
	} else {
		stop("target must be a character, logical, or numeric vector")
	}

	if (all(r == 0)) warning("No target states detected in 'target'")

	# bidirectional strength matrix
	B <- A + t(A)
	diag(B) <- 0

	# connectivity to targets and non-targets: sum_j B_ij * r_j or (1 - r_j)
	raw_to_target <- as.numeric(B %*% r)
	raw_to_nontarget <- as.numeric(B %*% (1 - r))

	if (normalize) {
		row_tot <- rowSums(B)
		score <- ifelse(row_tot > 0, raw_to_target / row_tot, 0)
		score_nt <- ifelse(row_tot > 0, raw_to_nontarget / row_tot, 0)
	} else {
		score <- raw_to_target
		score_nt <- raw_to_nontarget
	}

	state_names <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(n))
	is_target <- as.logical(r)

	out <- data.frame(
		state = state_names,
		score = score,
		score_nt = score_nt,
		raw_score = raw_to_target,
		raw_score_nt = raw_to_nontarget,
		is_target = is_target,
		stringsAsFactors = FALSE
	)

	# order: targets first (by proximity to non-targets, increasing),
	# then non-targets (by proximity to targets, decreasing)
	idx_target <- which(is_target)
	idx_nontarget <- which(!is_target)

	if (length(idx_target) > 0) {
		order_target <- idx_target[order(score_nt[idx_target], decreasing = FALSE)]
	} else {
		order_target <- integer(0)
	}
	if (length(idx_nontarget) > 0) {
		order_nontarget <- idx_nontarget[order(score[idx_nontarget], decreasing = TRUE)]
	} else {
		order_nontarget <- integer(0)
	}

	ord <- c(order_target, order_nontarget)
	out <- out[ord, , drop = FALSE]
	row.names(out) <- NULL
	out
}


# Compute average z-score per cluster for each feature
# feature_mat: matrix-like with rows = features, cols = cells
# cluster_dict: named vector mapping cell -> cluster
cluster_feature_zmeans <- function(feature_mat, cluster_dict) {
    if (is.null(colnames(feature_mat))) stop("feature_mat must have column names (cell IDs)")
    if (is.null(names(cluster_dict))) stop("cluster_dict must be a named vector mapping cell -> cluster")

    cells <- intersect(colnames(feature_mat), names(cluster_dict))
    if (length(cells) == 0) stop("No overlapping cells between feature_mat and cluster_dict")

    mat <- feature_mat[, cells, drop = FALSE]
    cl <- as.character(cluster_dict[cells])

    # per-feature z-score (center + scale across cells) using base `scale`
    # scale() works column-wise, so transpose twice to operate on rows
    z <- t(scale(t(mat), center = TRUE, scale = TRUE))
    # replace NaN produced by zero-variance rows with 0 so cluster means are 0
    z[is.na(z)] <- 0

    clusters <- unique(cl)
    out <- matrix(NA_real_, nrow = nrow(z), ncol = length(clusters), dimnames = list(rownames(z), clusters))
    for (k in clusters) {
        cols <- which(cl == k)
        if (length(cols) == 0) next
        out[, k] <- rowMeans(z[, cols, drop = FALSE], na.rm = TRUE)
    }
    as.data.frame(out, stringsAsFactors = FALSE)
}

get_marker_dat <- function(
	df_markers,
	state_order,
	cluster_dict,
	feature_mat,
	n_per_cluster = 10,
	min_log2fc = 0.25,
	max_q = 0.01,
	min_auc = 0.7
) {
	# 1. filter by single-cell DE + AUC
	df_filt <- df_markers %>%
		filter(cluster %in% state_order) %>%
		filter(p_val_adj < max_q) %>%
		filter(avg_log2FC > min_log2fc)

	genes_use <- unique(df_filt$gene)

	# 2. pseudobulk z (genes x clusters)
	z_mat <- feature_mat[genes_use, , drop = FALSE] %>%
		cluster_feature_zmeans(cluster_dict)
	z_mat <- z_mat[, state_order, drop = FALSE]

	# long format for joining
	z_long <- z_mat %>%
		as.data.frame() %>%
		tibble::rownames_to_column("gene") %>%
		tidyr::pivot_longer(
			cols = -gene,
			names_to = "cluster",
			values_to = "z_pb"
		)

	# 3. for each gene, keep only the cluster with max pseudobulk z
	z_max <- z_long %>%
		group_by(gene) %>%
		filter(z_pb == max(z_pb)) %>%
		slice_head(n = 1) %>%   # in case of ties, keep one
		ungroup()

	# 4. intersect with DE table and pick top N per cluster
	top_markers <- df_filt %>%
		inner_join(z_max, by = c("gene", "cluster")) %>%
		mutate(cluster = factor(cluster, levels = state_order)) %>%
		group_by(cluster) %>%
		arrange(desc(avg_log2FC)) %>%
		slice_head(n = n_per_cluster) %>%
		ungroup()

	# 5. build df_z for plotting
	feature_order <- top_markers %>% 
         arrange(cluster, desc(avg_log2FC)) %>%
         pull(gene) %>% unique()

	z_sub <- z_mat[feature_order, state_order, drop = FALSE]

	df_z <- z_sub %>%
		tibble::rownames_to_column("feature") %>%
		reshape2::melt(id.var = "feature", variable.name = "cluster", value.name = "score") %>%
		mutate(feature = factor(feature, rev(unique(feature))))

	list(df_z = df_z, top_markers = top_markers)
}

