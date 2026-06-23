#' Drop root-adjacent singleton tips from a phylo tree.
#'
#' @param phy A rooted `phylo` object (ape) with `edge` and `tip.label` slots.
#'
#' @return A `phylo` with any tips directly attached to the root removed.
#'
#' @details Identifies the unique root as the node never appearing as a child,
#' then prunes tips that connect directly to it via `ape::drop.tip`.
drop_singletons <- function(phy) {
	# Remove tips that are directly connected to the root (singletons)
	# phy: an object of class 'phylo' from the ape package

  if (!inherits(phy, "phylo")) {
    stop("Input must be a phylo object.")
  }
  if (is.null(phy$edge) || is.null(phy$tip.label)) {
    stop("Malformed phylo object.")
  }

  # Find the root node (the node that is not a child in any edge)
  all_nodes <- unique(as.vector(phy$edge))
  child_nodes <- unique(phy$edge[,2])
  root_node <- setdiff(all_nodes, child_nodes)
  if (length(root_node) != 1) {
    stop("Could not uniquely identify root node.")
  }
  root_node <- root_node[1]

  # Find all tips (leaves) directly connected to the root
  tip_nums <- seq_along(phy$tip.label)
  tips_connected_to_root <- phy$edge[phy$edge[,1] == root_node & phy$edge[,2] %in% tip_nums, 2]

  # Remove these tips
  if (length(tips_connected_to_root) > 0) {
    phy <- ape::drop.tip(phy, phy$tip.label[tips_connected_to_root])
  }
  phy
}


#' Keep a safe subset of tips present in the tree.
#'
#' @param phy A `phylo` object.
#' @param tips Character vector of tip labels to retain; labels not present are ignored.
#'
#' @return A `phylo` pruned to the intersection of `tips` and existing labels.
keep_tip_safe = function(phy, tips) {
    tips = intersect(tips, phy$tip.label)
    phy_sub = ape::keep.tip(phy, tips)
    return(phy_sub)
}



#' Compute Internal Node Scores
#'
#' Calculates scores for internal nodes based on the average scores of their descendant tips.
#' It computes a z-score and p-value for each internal node by comparing the clade mean
#' to the global mean and standard deviation of tip scores, accounting for clade size.
#'
#' @param tree A `phylo` object. Should have `tip.label` and optionally `node.label`.
#' @param tip_scores A named numeric vector of scores for the tips. Names must match `tree$tip.label`.
#' @param singletons Logical. If `TRUE`, includes rows for tip nodes (singletons) in the output,
#'   inheriting scores from their parent internal node. Default is `FALSE`.
#' @param min_clade_size Integer. Minimum clade size to include in the multiple-testing
#'   correction. Nodes smaller than this contribute raw p/z but are excluded from the
#'   BH denominator. Default is 1.
#'
#' @return A data frame with columns:
#'   \item{node}{Node ID (integer).}
#'   \item{label}{Node label (character).}
#'   \item{clade_size}{Number of tips in the clade.}
#'   \item{clade_mean}{Mean score of tips in the clade.}
#'   \item{avg_delta}{Difference between clade mean and global mean.}
#'   \item{z_value}{Z-score of the node score.}
#'   \item{z_adj}{Z-score adjusted for multiple testing (FDR).}
#'   \item{p_value}{Two-sided p-value based on the z-score.}
#'
#' @export
compute_internal_node_scores <- function(tree, tip_scores, singletons = FALSE, min_clade_size = 1) {
	
	# align scores to tree tip order
	if (is.null(names(tip_scores))) {
		stop("tip_scores must be a *named* vector with names = tree$tip.label")
	}
	tip_scores <- tip_scores[tree$tip.label]
	if (any(is.na(tip_scores))) {
		stop("Missing scores for some tips after aligning to tree$tip.label")
	}
	
	n_tip <- length(tree$tip.label)
	int_nodes <- (n_tip + 2):(n_tip + tree$Nnode)
	
	# clade sizes (tips only) using TreeTools
	clade_sizes <- TreeTools::CladeSizes(
		tree,
		internal = FALSE,           # count tips only
		nodes    = int_nodes
	)
	
	# descendants of each internal node (tip indices)
	desc_list <- phangorn::Descendants(
		tree,
		int_nodes,
		type = "tips"
	)
	
	# mean score within each clade
	clade_means <- vapply(
		desc_list,
		function(tips) mean(tip_scores[tips], na.rm = TRUE),
		numeric(1)
	)
	
	# global mean / sd (over all tips)
	global_mean <- mean(tip_scores, na.rm = TRUE)
	global_sd   <- stats::sd(tip_scores, na.rm = TRUE)
	
	# deviation from global mean
	avg_delta <- clade_means - global_mean
	
	# size-aware SE and z
	se <- global_sd / sqrt(clade_sizes)
	z  <- avg_delta / se
	
	# two-sided p if you want it
	p  <- 2 * stats::pnorm(-abs(z))

    # adjust for multiple testing (only nodes at/above size threshold)
    test_mask <- clade_sizes >= min_clade_size
    p_adj <- rep(NA_real_, length(p))
    if (any(test_mask)) {
        p_adj[test_mask] <- stats::p.adjust(p[test_mask], method = "BH")
    }

    z_adj <- rep(NA_real_, length(z))
    if (any(test_mask)) {
        z_adj[test_mask] <- sign(z[test_mask]) * stats::qnorm(p_adj[test_mask] / 2, lower.tail = FALSE)
    }
	
	out <- data.frame(
		node        = int_nodes,
		label       = if (!is.null(tree$node.label)) tree$node.label[int_nodes - n_tip] else NA_character_,
		clade_size  = as.integer(clade_sizes),
		clade_mean  = clade_means,
		avg_delta  = avg_delta,
		z_value     = z,
        z_adj       = z_adj,
		p_value     = p,
		row.names   = NULL,
		check.names = FALSE
	)

    if (singletons) {
        edge <- tree$edge
        parents <- edge[match(1:n_tip, edge[,2]), 1]
        idx <- match(parents, out$node)
        
        out_tips <- data.frame(
            node        = 1:n_tip,
            label       = tree$tip.label,
            clade_size  = 1L,
            clade_mean  = tip_scores,
            avg_delta  = out$avg_delta[idx],
            z_value     = out$z_value[idx],
            z_adj       = out$z_adj[idx],
            p_value     = out$p_value[idx],
            row.names   = NULL,
            check.names = FALSE
        )
        out <- rbind(out, out_tips)
    }
	
	out
}

#' Propagate node scores down the tree and return a named vector
#'
#' @param tree       ape::phylo
#' @param score_dict named numeric vector of scores (per node or per label)
#' @param threshold  absolute cutoff; abs(score) >= threshold starts/overrides a clade
#'
#' @return named numeric vector of propagated scores (names = tip + node labels)
propagate_scores_down_dict <- function(tree, score_dict, threshold) {
	if (is.null(names(score_dict))) {
		stop("score_dict must be a named numeric vector.")
	}

	n_tip  <- length(tree$tip.label)
	n_node <- tree$Nnode
	N      <- n_tip + n_node

	# all labels (tips then internals)
	all_labels <- c(
		tree$tip.label,
		if (!is.null(tree$node.label)) tree$node.label else as.character((n_tip + 1):N)
	)

	# map score_dict to node indices 1..N
	label_to_id <- setNames(seq_len(N), all_labels)
	names_sd    <- names(score_dict)

	if (all(names_sd %in% as.character(seq_len(N)))) {
		# treat names as node IDs
		node_ids <- as.integer(names_sd)
	} else {
		dup_labels <- unique(all_labels[duplicated(all_labels)])
		if (length(dup_labels) > 0) {
			ambig <- intersect(names_sd, dup_labels)
			if (length(ambig) > 0) {
				stop("Ambiguous label mapping: tree has duplicated tip/node labels")
			}
		}

		# treat names as labels
		if (!all(names_sd %in% names(label_to_id))) {
			miss <- setdiff(names_sd, names(label_to_id))
			stop("Some names in score_dict do not match any tip/node label: ",
				paste(head(miss, 10), collapse = ", "))
		}
		node_ids <- as.integer(label_to_id[names_sd])
	}

	score_node <- rep(NA_real_, N)
	score_node[node_ids] <- as.numeric(score_dict)

    pass_thr <- !is.na(score_node) & abs(score_node) >= threshold

	# parent/children structure
	E <- tree$edge
	parent <- rep(NA_integer_, N)
	parent[E[, 2]] <- E[, 1]
	children <- split(E[, 2], E[, 1])
	names(children) <- as.character(names(children))

	root <- setdiff(E[, 1], E[, 2])[1]

	# top-down propagation
	score_prop <- rep(NA_real_, N)

	walk <- function(node, current_score) {
		if (pass_thr[node]) {
			new_score <- score_node[node]
		} else {
			new_score <- current_score
		}
		score_prop[node] <<- new_score

		ch <- children[[as.character(node)]]
		if (!length(ch)) return(invisible(NULL))
		for (v in ch) {
			walk(v, new_score)
		}
	}

	walk(root, NA_real_)
	setNames(score_prop, all_labels)
}

#' Propagate scores but keep the maximum encountered along the root path
#'
#' @inheritParams propagate_scores_down_dict
#' @return named numeric vector where each node/tip stores the largest score
#'   (by raw value) observed on the path from that node to the root among
#'   threshold-passing nodes.
propagate_scores_down_max <- function(tree, score_dict, threshold = 0) {
	if (is.null(names(score_dict))) {
		stop("score_dict must be a named numeric vector.")
	}

	n_tip  <- length(tree$tip.label)
	n_node <- tree$Nnode
	N      <- n_tip + n_node

	all_labels <- c(
		tree$tip.label,
		if (!is.null(tree$node.label)) tree$node.label else as.character((n_tip + 1):N)
	)

	label_to_id <- setNames(seq_len(N), all_labels)
	names_sd    <- names(score_dict)

	if (all(names_sd %in% as.character(seq_len(N)))) {
		node_ids <- as.integer(names_sd)
	} else {
		dup_labels <- unique(all_labels[duplicated(all_labels)])
		if (length(dup_labels) > 0) {
			ambig <- intersect(names_sd, dup_labels)
			if (length(ambig) > 0) {
				stop("Ambiguous label mapping: tree has duplicated tip/node labels")
			}
		}

		if (!all(names_sd %in% names(label_to_id))) {
			miss <- setdiff(names_sd, names(label_to_id))
			stop("Some names in score_dict do not match any tip/node label: ",
				paste(head(miss, 10), collapse = ", "))
		}
		node_ids <- as.integer(label_to_id[names_sd])
	}

	score_node <- rep(NA_real_, N)
	score_node[node_ids] <- as.numeric(score_dict)

	pass_thr <- !is.na(score_node) & abs(score_node) >= threshold

	E <- tree$edge
	children <- split(E[, 2], E[, 1])
	names(children) <- as.character(names(children))
	root <- setdiff(E[, 1], E[, 2])[1]

	score_prop <- rep(NA_real_, N)

	walk <- function(node, best_score) {
		candidate <- best_score
		if (pass_thr[node]) {
			if (is.na(candidate) || abs(score_node[node]) > abs(candidate)) {
				candidate <- score_node[node]
			}
		}
		score_prop[node] <<- candidate

		ch <- children[[as.character(node)]]
		if (!length(ch)) return(invisible(NULL))
		for (v in ch) {
			walk(v, candidate)
		}
	}

	walk(root, NA_real_)
	setNames(score_prop, all_labels)
}


#' Majority element (mode)
#'
#' Return the most frequent value in a vector. In the case of ties,
#' the first value to appear in `x` (i.e., the order preserved by
#' `unique(x)`) is chosen.
#'
#' @param x A vector (numeric, character, factor, or logical).
#' @return A single value: the most frequent entry in `x`.
#' @examples
#' majority(c(1, 2, 2, 3))
#' majority(c("A", "B", "A"))
#' majority(factor(c("A", "B", "A")))
#' @export
majority <- function(x) {
  uniq_v <- unique(x)
  output <- uniq_v[which.max(tabulate(match(x, uniq_v)))]
  return(output)
}


# Generalized KNN neighbor names from an arbitrary-dimensional embedding.
# Input:
#   embed         : matrix/data.frame of size n x d (d arbitrary) with rownames = cell IDs
#   k             : number of neighbors to return per cell
#   include_self  : if TRUE, include the cell itself among its neighbors
# Output:
#   character matrix of size n x k_eff (rows = cells, columns = "nn1", "nn2", ...)
knn_embed_neighbors <- function(embed, k = 15, include_self = FALSE) {
  # normalize input
  if (is.data.frame(embed)) embed <- as.matrix(embed)
  if (!is.matrix(embed)) stop("`embed` must be a matrix or data.frame.")
  if (is.null(rownames(embed))) stop("`embed` must have rownames (cell IDs).")
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("Please install the 'RANN' package.")

  n <- nrow(embed)
  if (n == 0L) {
    out <- matrix(character(0), 0, 0)
    return(out)
  }
  if (n == 1L) {
    if (include_self && k >= 1) {
      out <- matrix(rownames(embed), 1, 1)
      rownames(out) <- rownames(embed)
      colnames(out) <- "nn1"
      return(out)
    } else {
      out <- matrix(character(0), 1, 0)
      rownames(out) <- rownames(embed)
      return(out)
    }
  }

  # Cap k to available neighbors
  k_eff <- max(0L, min(as.integer(k), n - as.integer(!include_self)))
  k_query <- k_eff + as.integer(!include_self)  # add 1 to drop self if needed
  if (k_eff == 0L) {
    out <- matrix(character(0), n, 0)
    rownames(out) <- rownames(embed)
    return(out)
  }

  nn_idx <- RANN::nn2(data = embed, query = embed, k = k_query)$nn.idx
  if (!include_self) nn_idx <- nn_idx[, -1, drop = FALSE]  # drop self

  res <- matrix(rownames(embed)[nn_idx],
                nrow = n,
                ncol = ncol(nn_idx))
  rownames(res) <- rownames(embed)
  colnames(res) <- paste0("nn", seq_len(ncol(res)))
  res
}


compute_new_wnn = function(seu, cluster_res = 1) {
    # 1) Preprocess RNA with standard normalization + PCA
    DefaultAssay(seu) <- "RNA"
    seu = SCTransform(seu, verbose = FALSE)
    DefaultAssay(seu) <- "SCT"
    seu <- RunPCA(seu, reduction.name = 'pca.RNA', features = VariableFeatures(seu), verbose = FALSE)

    DefaultAssay(seu) <- "ATAC"
    seu <- RunTFIDF(seu, verbose = FALSE)
    seu <- FindTopFeatures(seu, min.cutoff = 'q0', verbose = FALSE)
    seu <- RunSVD(seu, reduction.name = 'lsi.ATAC', verbose = FALSE)
    
    # 2) Build the WNN graph on RNA PCA and ATAC LSI reductions
    seu <- FindMultiModalNeighbors(
      seu,
      reduction.list = list(
        RNA     = "pca.RNA",
        ATAC = "lsi.ATAC"
      ),
      dims.list = list(
        RNA     = 1:30, 
        ATAC = 2:30
      ),
      verbose = FALSE
    )

    # 3) Cluster cells based on the WNN graph
    seu <- FindClusters(seu, graph.name = "wsnn", algorithm = 1, 
        resolution = cluster_res, verbose = FALSE,
        cluster.name = paste0("wsnn_atac_rna_res.", cluster_res))
    
    # 4) Run UMAP on the weighted-NN graph for visualization
    seu <- RunUMAP(
      object    = seu,
      nn.name   = "weighted.nn",               # the WNN graph stored by FindMultiModalNeighbors
      reduction = "weighted.nn",               # name of the WNN graph
      reduction.name = "wnn.umap.RNA_ATAC",     # your custom UMAP slot
      verbose = FALSE
    )

    return(seu)
}


compute_comdist <- function(
    tree,
    cluster_dict,
    B = 0,                      # number of bootstrap replicates (0 = no CI, return matrix as before)
    level = 0.95,               # CI level
    abundance.weighted = FALSE,  # pass-through to picante::comdist
    ncores = 1,                 # POSIX only for parallel::mclapply
    seed = 1L,                  # RNG seed base
    keep_boot = FALSE           # optionally return the 3D bootstrap array
) {

    RhpcBLASctl::blas_set_num_threads(ncores)

    # Pairwise phylogenetic distances among tips (unit branch lengths, unweighted topology)
    tree$edge.length <- rep(1, nrow(tree$edge))
    M <- ape::cophenetic.phylo(tree)
    print(dim(M))

    # subset to tips present in cluster_dict
    tips = intersect(tree$tip.label, names(cluster_dict))
    M = M[tips, tips]

    # Community-by-tip matrix X (rows = clusters/communities, cols = tips/cells)
    X <- PATH::catMat(cluster_dict[tips]) |> t() |> as.matrix()
    colnames(X) <- tips

    # Observed COMDIST matrix (as a numeric matrix for easy CI work)
    obs <- as.matrix(picante::comdist(X, M, abundance.weighted = abundance.weighted))

    # If no bootstrap requested, return long-form data.frame with obs and NA CIs
    if (B <= 0) {
        df <- reshape2::melt(obs) %>%
            setNames(c("com_x", "com_y", "obs")) %>%
            dplyr::mutate(lower = NA_real_, upper = NA_real_)
        return(df)
    }

    # --- Bootstrap helper: resample each community row with replacement via multinomial ---
    resample_X <- function(Xrow) {
        s <- sum(Xrow)
        if (s <= 0) return(Xrow)                 # empty community → stay empty
        p <- Xrow / s
        as.vector(stats::rmultinom(n = 1, size = s, prob = p))
    }

    # Run B replicates
    mats <- lapply(seq_len(B), function(b) {
        # (Optionally) make RNG per-replicate reproducible
        if (!is.null(seed)) set.seed(seed + b)
        Xb <- t(apply(X, 1L, resample_X))
        rownames(Xb) <- rownames(X)
        colnames(Xb) <- colnames(X)
        as.matrix(picante::comdist(Xb, M, abundance.weighted = abundance.weighted))
    })

    # Stack into an array: [community, community, bootstrap]
    K <- nrow(obs)
    labs <- rownames(obs)
    arr <- array(NA_real_, dim = c(K, K, B), dimnames = list(labs, labs, NULL))
    for (b in seq_len(B)) arr[,,b] <- mats[[b]]

    # Pairwise bootstrap CIs
    alpha <- (1 - level) / 2
    lower <- apply(arr, c(1, 2), stats::quantile, probs = alpha, na.rm = TRUE)
    upper <- apply(arr, c(1, 2), stats::quantile, probs = 1 - alpha, na.rm = TRUE)

    # Build long-form data.frame with obs, lower, upper per community pair
    df_obs   <- reshape2::melt(obs)   %>% setNames(c("com_x", "com_y", "obs"))
    df_lower <- reshape2::melt(lower) %>% setNames(c("com_x", "com_y", "lower"))
    df_upper <- reshape2::melt(upper) %>% setNames(c("com_x", "com_y", "upper"))

    df <- df_obs %>%
        dplyr::left_join(df_lower, by = c("com_x", "com_y")) %>%
        dplyr::left_join(df_upper, by = c("com_x", "com_y"))

    return(df)
}

#' Add "Node<n>" labels to the internal nodes of a phylo tree
#'
#' @param tree A phylo object
#' @param prefix Character prefix for node names (default "Node")
#' @return The same phylo object, with tree$node.label set to prefix + node numbers
add_node_names <- function(tree, prefix = "Node", start_from_tip = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylo object")
  }
  ntip  <- length(tree$tip.label)
  nnode <- tree$Nnode

  # internal node IDs run from ntip+1 to ntip+nnode
  if (start_from_tip) {
    ids <- seq(ntip + 1, ntip + nnode)
  } else {
    ids <- seq(1, nnode)
  }

  # assign labels
  tree$node.label <- paste0(prefix, ids)

  return(tree)
}

#' Compute root-to-node depths for internal nodes
#'
#' Returns the root-to-node depth (height) for every internal node in `tree`,
#' computed from branch lengths via `ape::node.depth.edgelength()`.
#'
#' @param tree A `phylo` object with branch lengths.
#' @param include_root If FALSE, omit the root node from the returned depths.
#' @return Named numeric vector of internal-node depths (names are node IDs).
get_node_depths <- function(tree, include_root = TRUE) {
    if (is.null(tree$edge.length)) stop("tree has no branch lengths (edge.length is NULL).")

	depth <- ape::node.depth.edgelength(tree)
	n_tip <- length(tree$tip.label)
	internal_ids <- n_tip + seq_len(tree$Nnode)
    if (!include_root) {
        root <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
        internal_ids <- internal_ids[internal_ids != root]
    }
	setNames(depth[internal_ids], internal_ids)
}

#' Compute internal node ages relative to the tips
#'
#' Defines node age as the difference between the mean tip height and the node
#' height, i.e. $\bar{h}_{tips} - h_{node}$ where heights are root-to-node
#' distances computed from branch lengths.
#'
#' @param tree A `phylo` object with branch lengths.
#' @return Named numeric vector of ages for internal nodes (names are node IDs).
compute_internal_node_ages <- function(tree) {
    if (is.null(tree$edge.length)) stop("tree has no branch lengths (edge.length is NULL).")

    depth <- ape::node.depth.edgelength(tree)
    n_tip <- length(tree$tip.label)
    internal_ids <- n_tip + seq_len(tree$Nnode)

    mean_tip_height <- mean(depth[seq_len(n_tip)])
    ages <- mean_tip_height - depth[internal_ids]
    setNames(ages, internal_ids)
}

#' Compute Shannon diversity (entropy) for a categorical vector
#'
#' Computes Shannon entropy of the empirical distribution of values in `vec`:
#' \deqn{H = - \sum_i p_i \log(p_i)}
#' and converts it to the requested log base via \eqn{H / \log(base)}.
#'
#' @param vec Vector of categories (character, factor, integer, etc.). `NA` values
#'   are ignored (as per `table()` default behavior).
#' @param base Log base for the returned entropy. Use `exp(1)` for nats
#'   (default), `2` for bits, or `10` for bans.
#'
#' @return A single numeric value: Shannon entropy of `vec` in the chosen base.
calc_sdi <- function(vec, base = exp(1)) {
	tbl <- table(vec)
	p <- as.numeric(tbl) / sum(tbl)
	H_ln <- -sum(p * log(p))              # natural log
	H <- H_ln / log(base)                 # convert to chosen base (e, 2, 10)
	return(H)
}

#' Build a named Brewer palette with optional cycling behavior
#'
#' Creates up to `n` visually distinct colors by sampling (and, when needed,
#' interpolating) from an RColorBrewer palette, optionally reordering each
#' block in a zig-zag pattern and progressively lightening subsequent cycles.
#'
#' @param n Integer number of colors to generate.
#' @param labels Optional character vector of names to assign to the colors. If
#'   `NULL`, sequential integers are used.
#' @param pal RColorBrewer palette name provided to `RColorBrewer::brewer.pal`.
#' @param alternating Logical; when `TRUE`, reorders each chunk via a
#'   zig-zag pattern to maximize separation of adjacent hues.
#' @param cycle_len Positive integer forcing palette generation in repeating
#'   chunks of this size; `0` (default) disables chunking.
#' @param cycle_shift Scalar in `[0, 1]` controlling how strongly later chunks
#'   blend toward white to maintain distinguishability.
#'
#' @return Named character vector of hex colors.
make_clade_pal <- function(
	n,
	labels = NULL,
	pal = "Set3",
	alternating = FALSE,
	cycle_len = 0L,
	cycle_shift = 0.15
) {

	max_n <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
	base_cols <- RColorBrewer::brewer.pal(max_n, pal)

	zigzag_idx <- function(len) {
		if (len <= 2L) return(seq_len(len))
		ord <- as.vector(rbind(seq_len(len), rev(seq_len(len))))
		ord <- ord[ord <= len]
		ord[!duplicated(ord)]
	}

	blend_with_white <- function(cols, frac) {
		if (frac <= 0) return(cols)
		frac <- pmin(frac, 0.9)
		rgb_mat <- grDevices::col2rgb(cols) / 255
		mix <- (1 - frac) * rgb_mat + frac
		grDevices::rgb(mix[1, ], mix[2, ], mix[3, ])
	}

	if (cycle_len <= 0L) {
		cols <- if (n <= max_n) {
			base_cols[seq_len(n)]
		} else {
			grDevices::colorRampPalette(base_cols)(n)
		}
		if (alternating) {
			cols <- cols[zigzag_idx(length(cols))]
		}
	} else {
		chunk_size <- max(1L, min(cycle_len, max_n))
		base_chunk <- if (chunk_size <= max_n) {
			base_cols[seq_len(chunk_size)]
		} else {
			grDevices::colorRampPalette(base_cols)(chunk_size)
		}

		cols <- character(0)
		remaining <- n
		cycle_idx <- 1L
		while (remaining > 0L) {
			take <- min(chunk_size, remaining)
			block <- base_chunk
			if (take < length(block)) {
				block <- block[seq_len(take)]
			}
			if (alternating && take > 1L) {
				block <- block[zigzag_idx(take)]
			}
			if (cycle_idx > 1L) {
				frac <- min(0.85, (cycle_idx - 1L) * cycle_shift)
				block <- blend_with_white(block, frac)
			}
			cols <- c(cols, block)
			remaining <- remaining - take
			cycle_idx <- cycle_idx + 1L
		}
	}

	if (is.null(labels)) {
		labels <- as.character(seq_len(n))
	} else {
		labels <- as.character(labels)
	}

	names(cols) <- labels
	cols
}

#' Generate a color palette with distinct colors
make_named_polychrome_pal <- function(n, labels = NULL, seed = 1L) {

	if (is.null(labels)) {
		labels <- as.character(seq_len(n))
	} else {
		labels <- as.character(labels)
		if (length(labels) != n) {
			stop("length(labels) must equal n when labels are provided")
		}
	}

	default_seeds <- c(
		"#4E79A7",  # muted blue
		"#59A14F",  # muted green
		"#E15759",  # muted red
		"#F28E2B",  # muted orange
		"#B07AA1",  # muted purple
		"#76B7B2"   # muted teal
	)

	set.seed(as.integer(seed))

	cols <- Polychrome::createPalette(n, seedcolors = default_seeds)
	names(cols) <- labels
	cols
}


#' k-NN category counts on a phylogeny
#'
#' Given a phylo object and a named vector mapping tip labels to categories,
#' compute, for each tip, the counts (or proportions) of categories among its
#' k nearest neighbors by patristic distance.
#'
#' @param phy		An `ape::phylo` object.
#' @param cat_dict	Named vector (character or factor). Names must be tip labels; values are categories.
#' @param k		Integer; number of nearest neighbors to use (default 20).
#' @param include_self Logical; if TRUE, a tip may be its own neighbor (default FALSE).
#' @param cat_exclude Optional vector of category labels to *exclude* from the neighbor pool.
#' @param return_prop Logical; if TRUE, return proportions instead of counts (default FALSE).
#' @param ncores	Integer number of cores for mclapply (POSIX only).
#' @param seed		Integer RNG seed base for tie-breaking (per-tip uses seed+i for stability).
#' @param ties		"random" (sample within equal-distance ties to keep exactly k) or
#' 				"include" (include whole tie groups; may exceed k).
#' @param return_neighbors Logical; if TRUE, appends a list-column `neighbors` with the
#' 				character vector of selected neighbor tip names for each focal tip (default FALSE).
#' @return A data.frame with one row per tip and one column per category; when
#' 	`return_neighbors = TRUE`, also contains a list-column `neighbors`.
#' @examples
#' # counts <- knn_category_counts(phy, cat_dict, k = 15)
#' # props  <- knn_category_counts(phy, cat_dict, k = 15, return_prop = TRUE)
#' # with neighbors:
#' # res <- knn_category_counts(phy, cat_dict, k = 15, return_neighbors = TRUE)
knn_category_counts <- function(phy, cat_dict, k = 20, include_self = FALSE, cat_exclude = NULL,
	return_prop = TRUE, ncores = 1, seed = 0,
	ties = c("random", "include"), return_neighbors = FALSE) {
		ties <- match.arg(ties)
		if (!inherits(phy, "phylo")) stop("`phy` must be a phylo object")
		if (is.null(names(cat_dict))) stop("`cat_dict` must be a named vector with tip labels as names")
		# Align category vector to tree tips
		tips <- phy$tip.label
		cats_all <- cat_dict[tips]

		# Determine category universe (columns)
		cat_levels <- sort(unique(cats_all))
		
		# Distances (patristic)
		D <- NULL
		phy$edge.length = rep(1, nrow(phy$edge))  # unweighted
		D <- ape::cophenetic.phylo(phy)
		D <- D[tips, tips, drop = FALSE]

		if (!is.null(cat_exclude)) {
			exclude_idx <- which(cats_all %in% cat_exclude)
			if (length(exclude_idx) > 0) {
				D[ , exclude_idx] <- Inf
			}
		}

		# Parallel over tips: for each tip, compute k-NN category counts or proportions
		res_list <- parallel::mclapply(seq_along(tips), mc.cores = ncores, function(i) {
			set.seed(seed+i) # per-tip seed for tie-breaking; this is important otherwise it oversamples certain tips
			# initialize zero vector for this row
			row <- setNames(numeric(length(cat_levels)), as.character(cat_levels))
			di <- D[i, ]
			names(di) <- tips
			self_name <- if (include_self) NULL else tips[i]

			if (ties == "include") {
				sel_names <- knn_include_ties(di, k = k, self = self_name) 
			} else {
				sel_names <- knn_random_ties(di, k = k, self = self_name)
			}

			if (length(sel_names)) {
				nn_idx <- match(sel_names, tips)
				ci <- cats_all[nn_idx]
				if (length(ci)) {
					ct <- table(ci)
					row[names(ct)] <- as.numeric(ct)
				}
			}

			if (return_prop) {
				denom <- length(sel_names)
				if (denom > 0) row <- row / denom
			}
			list(row = row, neighbors = sel_names)
		})

		rows <- lapply(res_list, `[[`, "row")
		neighbors_list <- lapply(res_list, `[[`, "neighbors")

		out <- do.call(rbind, rows)
		out <- as.data.frame(out, stringsAsFactors = FALSE)
		rownames(out) <- tips
		if (isTRUE(return_neighbors)) {
			out$neighbors <- I(neighbors_list)
		}
		out
}

knn_random_ties <- function(dist_vec, k, self = NULL) {
	# dist_vec: named numeric vector of distances from a focal tip to all tips (may include the focal tip)
	# self: optional name of the focal tip to exclude
	if (!is.null(self) && self %in% names(dist_vec)) {
		dist_vec[self] <- Inf
	}

	# Work through distances in increasing order; when multiple neighbors share the
	# same distance, sample uniformly at random among that tie group.
	uniq_d <- sort(unique(dist_vec[is.finite(dist_vec)]))
	selected <- character(0)
	for (d in uniq_d) {
		group <- names(dist_vec)[which(dist_vec == d)]
		if (!length(group)) next
		if (length(selected) + length(group) <= k) {
			selected <- c(selected, sample(group, length(group)))
		} else {
			need <- k - length(selected)
			if (need > 0) {
				selected <- c(selected, sample(group, need))
			}
			break
		}
	}
	return(selected)
}

# Like knn_random_ties, but when the next distance tier would exceed k,
# include the entire tie group (so the returned set may be slightly larger than k).
# This avoids arbitrarily sampling within ties and instead keeps all equidistant neighbors.
knn_include_ties <- function(dist_vec, k, self = NULL) {
    # dist_vec: named numeric vector of distances from a focal tip to all tips (may include the focal tip)
    # self: optional name of the focal tip to exclude
    if (!is.null(self) && self %in% names(dist_vec)) {
        dist_vec[self] <- Inf
    }

    # consider only finite distances, in increasing order
    uniq_d <- sort(unique(dist_vec[is.finite(dist_vec)]))
    selected <- character(0)
    for (d in uniq_d) {
        group <- names(dist_vec)[which(dist_vec == d)]
        # always include the whole tie group
        selected <- c(selected, group)
        # stop once we have reached or exceeded k after including a full group
        if (length(selected) >= k) break
    }
    unique(selected)
}


# KNN smoothing on arbitrary embedding coords for any value column
knn_smooth_embedding <- function(
	df,
	value_col,                       # e.g. "myeloid"
	coord_cols = c("umap_x","umap_y"),
	k = 30,
	include_self = FALSE,
	weighted = FALSE,                # Gaussian distance weighting
	sigma = NULL,                    # bandwidth; if NULL, per-point median kNN dist
	out_col = NULL,                  # default "<value_col>_smooth"
	na_rm = TRUE                     # ignore NA neighbors when averaging
) {
	# checks
	if (!all(coord_cols %in% names(df))) stop("Missing coord_cols in df.")
	if (!value_col %in% names(df)) stop("Missing value_col in df.")
	X <- as.matrix(df[, coord_cols, drop = FALSE])
	y <- as.numeric(df[[value_col]])
	n <- nrow(df)
	if (n < 2) stop("Need at least 2 rows.")

	# neighbors
	k_search <- min(n, k + (!include_self))
	nn <- RANN::nn2(X, query = X, k = k_search)
	idx <- nn$nn.idx
	dist <- nn$nn.dists
	if (!include_self) {
		idx <- idx[, -1, drop = FALSE]
		dist <- dist[, -1, drop = FALSE]
	}
	k_eff <- ncol(idx)  # may be < k if n small

	# gather neighbor values into matrix (n x k_eff)
	Ymat <- matrix(y[idx], nrow = n)

	# weights
	if (weighted) {
		if (is.null(sigma)) {
			sigma <- matrixStats::rowMedians(dist)
			sigma <- pmax(sigma, 1e-12)  # avoid zero bandwidth
		}
		S <- matrix(sigma^2, nrow = n, ncol = k_eff)
		W <- exp(-(dist^2) / (2 * S))
		if (na_rm) {
			W[is.na(Ymat)] <- 0
			Yfill <- ifelse(is.na(Ymat), 0, Ymat)
			den <- pmax(rowSums(W), 1e-12)
			smooth <- rowSums(W * Yfill) / den
		} else {
			# if any NA present, result NA for that row
			smooth <- rowSums(W * Ymat) / rowSums(W)
		}
	} else {
		if (na_rm) {
			den <- rowSums(!is.na(Ymat))
			den <- pmax(den, 1L)
			Yfill <- ifelse(is.na(Ymat), 0, Ymat)
			smooth <- rowSums(Yfill) / den
		} else {
			smooth <- rowMeans(Ymat)
		}
	}

	# attach
	if (is.null(out_col)) out_col <- paste0(value_col, "_smooth")
	df[[out_col]] <- smooth
	df
}


get_discrete_colors <- function(categories, palette = "Set3") {

    # Ensure categories is a character vector and remove NA values
    categories <- as.character(categories)
    categories <- categories[!is.na(categories)]
    cats <- unique(categories)

    # Try to use RColorBrewer if available
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        pal_info <- RColorBrewer::brewer.pal.info
        if (!palette %in% rownames(pal_info)) {
            stop("Palette '", palette, "' not found in RColorBrewer.")
        }
        max_pal <- pal_info[palette, "maxcolors"]
        cols <- if (length(cats) <= max_pal) {
            RColorBrewer::brewer.pal(length(cats), palette)
        } else {
            grDevices::colorRampPalette(RColorBrewer::brewer.pal(max_pal, palette))(length(cats))
        }
    } else {
        warning("RColorBrewer not installed, falling back to rainbow().")
        cols <- grDevices::rainbow(length(cats))
    }

    # Create named mapping and map back to input
    col_map <- setNames(cols, cats)
    mapped_cols <- col_map[categories]

    col_map
}

is_transition <- function(ref, alt) {
	ref <- toupper(ref)
	alt <- toupper(alt)

	(ref == "A" & alt == "G") |
	(ref == "G" & alt == "A") |
	(ref == "C" & alt == "T") |
	(ref == "T" & alt == "C")
}

#' Convert mutation matrices to long format
#'
#' @param amat Alternate allele count matrix (variants x cells)
#' @param dmat Coverage/depth matrix (variants x cells)
#' @return Data frame in long format with columns: variant, cell, a, d
mut_mat_to_long <- function(amat, dmat) {
    # Convert sparse matrices to regular matrices if needed
    if (inherits(amat, "sparseMatrix")) {
        amat <- as.matrix(amat)
    }
    if (inherits(dmat, "sparseMatrix")) {
        dmat <- as.matrix(dmat)
    }

    a_long <- as.data.frame(amat) %>%
        tibble::rownames_to_column('variant') %>%
        reshape2::melt(id.var = 'variant', variable.name = 'cell', value.name = 'a')

    d_long <- as.data.frame(dmat) %>%
        tibble::rownames_to_column('variant') %>%
        reshape2::melt(id.var = 'variant', variable.name = 'cell', value.name = 'd')

    mut_dat_long <- full_join(
        a_long,
        d_long,
        by = join_by(variant, cell)
    ) %>%
    mutate(cell = as.character(cell))

    return(mut_dat_long)
}

#' Filter mutation data in long format
#' 
#' Removes homozygous variants and filters by minimum mean alt reads
#' 
#' @param long_matrix Data frame with columns: cell, variant, a (alt count), d (depth)
#' @param min_mean_alt Minimum mean alternate reads across cells with variant
#' @param min_frac_cells Fraction of cells threshold for identifying homozygous variants
#' @return Filtered data frame in long format
filter_mut_long = function(long_matrix, min_mean_alt = 1, min_frac_cells = 0.7) {

    hom_vars = long_matrix %>% 
        group_by(variant) %>% 
        mutate(vaf = a/d) %>%
        summarise(gvaf = sum(a)/sum(d), mean_vaf = mean(vaf), frac_cells = mean(a[d!=0]>0)) %>%
        filter(frac_cells > min_frac_cells) %>%
        pull(variant)

    message(paste0("Removing homozygous variants: ", paste(hom_vars, collapse = ", ")))

    mut_stats = long_matrix %>% 
        group_by(variant) %>% 
        filter(sum(a)>0) %>%
        summarise(mean_alt_reads = mean(a[a>0]), max_alt_reads = max(a[a>0]))

    variants_select = mut_stats %>% filter(mean_alt_reads > min_mean_alt) %>% pull(variant)

    mut_dat_full = long_matrix %>% 
        select(cell, variant, a, d) %>%
        filter(!variant %in% hom_vars) %>%
        filter(variant %in% variants_select) %>%
        group_by(variant) %>%
        filter(sum(a>0)>0) %>%
        ungroup()

    return(mut_dat_full)
}

get_clade_order <- function(clade_annot, phy_trim, flip = TRUE) {

    p_tmp <- ggtree(phy_trim, ladderize = TRUE, right = flip)
    cell_order <- p_tmp$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

    clade_order <- clade_annot %>%
      group_by(clade) %>%
      summarize(first_pos = min(match(cell, cell_order)), .groups = "drop") %>%
      arrange(first_pos) %>%
      pull(clade)

    return(clade_order)
}


#' Compute rank plot data for confidence diagnostics
#'
#' @param phy A phylo object with node.label containing confidence values
#' @param tree_name Optional name for tree grouping variable (default: NULL)
#'
#' @return Data frame with columns: node_rank, confidence, and optionally tree
compute_conf_rank_data <- function(phy, tree_name = NULL) {

    conf_values <- as.numeric(phy$node.label)

    df <- data.frame(
        node_rank = rank(-conf_values),  # Higher conf = lower rank
        confidence = conf_values
    )

    if (!is.null(tree_name)) {
        df$tree <- tree_name
        df <- df[, c("tree", "node_rank", "confidence")]
    }

    return(df)
}


#' Compute CDF data for confidence diagnostics
#'
#' @param phy A phylo object with node.label containing confidence values
#' @param conf_breaks Numeric vector of confidence thresholds to evaluate
#' @param tree_name Optional name for tree grouping variable (default: NULL)
#'
#' @return Data frame with columns: conf_cutoff, n_remaining, and optionally tree
compute_conf_cdf_data <- function(phy, conf_breaks = seq(0, 1, 0.01), tree_name = NULL) {

    conf_values <- as.numeric(phy$node.label)

    df <- data.frame(
        conf_cutoff = conf_breaks,
        n_remaining = vapply(conf_breaks, function(cutoff) {
            sum(conf_values >= cutoff)
        }, numeric(1))
    )

    if (!is.null(tree_name)) {
        df$tree <- tree_name
        df <- df[, c("tree", "conf_cutoff", "n_remaining")]
    }

    return(df)
}


#' Plot clade confidence diagnostics for phylo tree(s)
#'
#' @param phy A phylo object OR a named list of phylo objects, each with node.label containing confidence values
#' @param title Character string for plot title (default: NULL)
#' @param conf_breaks Numeric vector of confidence thresholds to evaluate (default: seq(0, 1, 0.01))
#' @param colors Named vector of colors for each tree when phy is a list (default: auto-assigned)
#'
#' @return A ggplot object of the clade-retention curve (remaining internal nodes vs confidence cutoff).
#'
#' @details
#' The plot shows how many internal nodes (clades) remain after filtering at
#' different confidence thresholds, analogous to a survival curve. Multiple
#' trees can be overlaid with distinct colors for direct comparison.
#'
#' @examples
#' # Single tree
#' phy <- read.tree("tree_annotated.newick")
#' plot_clade_conf_diagnostic(phy, title = "Sample Tree")
#'
#' # Multiple trees for comparison
#' phy1 <- read.tree("tree1_annotated.newick")
#' phy2 <- read.tree("tree2_annotated.newick")
#' plot_clade_conf_diagnostic(list("Version 1" = phy1, "Version 2" = phy2), title = "Comparison")
plot_clade_conf_diagnostic <- function(phy, title = NULL, conf_breaks = seq(0, 1, 0.01), colors = NULL) {

    # Normalize input to list format
    is_single_tree <- inherits(phy, "phylo")

    if (is_single_tree) {
        phy_list <- list(tree = phy)
    } else {
        phy_list <- phy
        if (!is.list(phy_list) || length(phy_list) == 0) {
            stop("phy must be a phylo object or a non-empty list of phylo objects")
        }
        if (is.null(names(phy_list))) {
            names(phy_list) <- paste0("Tree", seq_along(phy_list))
        }
    }

    # Default colors if not provided (only used for multiple trees)
    if (is.null(colors)) {
        colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")[seq_along(phy_list)]
        names(colors) <- names(phy_list)
    }

	# Compute data for all trees (CDF only)
	all_cdf <- lapply(names(phy_list), function(tree_name) {
		compute_conf_cdf_data(phy_list[[tree_name]], conf_breaks, tree_name)
	})

	df_cdf <- do.call(rbind, all_cdf)

	# Build CDF plot (single ggplot call)
	p_cdf <- ggplot(df_cdf, aes(x = conf_cutoff, y = n_remaining)) +
        geom_line(aes(color = if (is_single_tree) NULL else tree), linewidth = 0.8) +
        theme_bw() +
        scale_x_continuous(breaks = seq(0, 1, 0.2)) +
        labs(
            x = "Confidence cutoff",
            y = "# of remaining clades",
            color = if (is_single_tree) NULL else "Tree",
            title = if (!is.null(title)) paste0(title, ": Clade retention") else "Clade retention"
        ) +
        theme(
            panel.grid.minor = element_blank(),
			legend.position = if (is_single_tree) "none" else "bottom"
        )

    # Add color scale and points for appropriate cases
    if (!is_single_tree) {
        p_cdf <- p_cdf + scale_color_manual(values = colors)
    } else {
        p_cdf <- p_cdf + geom_point(size = 0.5, alpha = 0.3)
    }

	if (!is_single_tree) {
		p_cdf <- p_cdf + theme(legend.position = "bottom")
	}

	return(p_cdf)
}

#' Test clade x category associations using Fisher's exact test
#'
#' @param cat_dict Named vector mapping cell names to categories (e.g., clusters)
#' @param clone_dict Named vector mapping cell names to clones/clades
#'
#' @return Data frame with Fisher test results including p-values, odds ratios, and confidence intervals
#'
#' @details
#' For each clone-category pair, performs a 2x2 Fisher's exact test comparing:
#' - Cells in clone with category vs cells in clone without category
#' - Cells not in clone with category vs cells not in clone without category
#'
#' Returns results sorted by adjusted p-value (Benjamini-Hochberg correction).
#'
#' @examples
#' clone_dict <- setNames(clone_annot$clade, clone_annot$cell)
#' cluster_dict <- setNames(seu$seurat_clusters, colnames(seu))
#' fisher_res <- test_feature_fisher_2x2(cluster_dict, clone_dict)
test_feature_fisher_2x2 <- function(cat_dict, clone_dict, alternative = "two.sided") {
    # Filter clone_dict to only include cells that are in cat_dict
    common = intersect(names(clone_dict), names(cat_dict))
    clone_dict = clone_dict[common]
    cat_dict = cat_dict[common]

    # Get unique values
    categories <- sort(unique(cat_dict))
    clones <- sort(unique(clone_dict))

    # Initialize results list
    results <- list()

    # Loop through each clone and category combination
    for (clone in clones) {
        for (category in categories) {
            # Get cells in current clone
            clone_cells <- names(clone_dict)[clone_dict == clone]
            all_clone_cells <- names(cat_dict) %in% clone_cells

            # Create 2x2 contingency table
            # Rows: in clone vs not in clone
            # Columns: has category vs doesn't have category
            in_clone_with_cat <- sum(cat_dict[all_clone_cells] == category)
            in_clone_without_cat <- sum(cat_dict[all_clone_cells] != category)
            not_in_clone_with_cat <- sum(cat_dict[!all_clone_cells] == category)
            not_in_clone_without_cat <- sum(cat_dict[!all_clone_cells] != category)

            # Create 2x2 table
            cont_table <- matrix(
                c(in_clone_with_cat, not_in_clone_with_cat,
                  in_clone_without_cat, not_in_clone_without_cat),
                nrow = 2,
                dimnames = list(
                    c("in_clone", "not_in_clone"),
                    c("has_category", "no_category")
                )
            )

            # Perform Fisher test
            fisher_result <- fisher.test(cont_table, alternative = alternative)

            # Store results
            results[[paste0(clone, "_", category)]] <- data.frame(
                clone = clone,
                category = category,
                in_clone_with_cat = in_clone_with_cat,
                in_clone_without_cat = in_clone_without_cat,
                not_in_clone_with_cat = not_in_clone_with_cat,
                not_in_clone_without_cat = not_in_clone_without_cat,
                p.value = fisher_result$p.value,
                odds_ratio = fisher_result$estimate,
                conf_lower = fisher_result$conf.int[1],
                conf_upper = fisher_result$conf.int[2],
                stringsAsFactors = FALSE
            )
        }
    }

    # Combine all results
    results_df <- do.call(rbind, results)
    rownames(results_df) <- NULL

    # Add adjusted p-values
    results_df$p.adj <- p.adjust(results_df$p.value, method = "BH")

    # Sort by adjusted p-value
    results_df <- results_df[order(results_df$p.adj, results_df$p.value), ]

    return(results_df)
}
