test_feature_regression_covar <- function(feature_x_cell, cell_values, df_covar, ncores = 1) {
    feature_x_cell = as.matrix(feature_x_cell)
	if (is.null(names(cell_values))) stop("cell_values must be a named numeric vector")
	cells <- Reduce(intersect, list(colnames(feature_x_cell), names(cell_values)))
    message(glue::glue('{length(cells)} cells left after intersecting'))

	all_features <- rownames(feature_x_cell)

	results_list <- mclapply(all_features, mc.cores = ncores, function(feature) {

        df <- data.frame(
            cell = cells,
            y = feature_x_cell[feature, cells, drop = TRUE],
            x = cell_values[cells]
        ) %>%
        cbind(df_covar[cells,,drop = FALSE])
        
        # build formula
        covar_names <- setdiff(colnames(df_covar), "cell")
        formula_str <- paste("y ~ x +", paste(covar_names, collapse = " + "))
        
		fit <- lm(as.formula(formula_str), data = df)
		
		# Extract summary statistics
		summary_fit <- summary(fit)
		coef_table <- summary_fit$coefficients
		
		# Extract beta coefficient and p-value for the 'x' variable
		beta <- coef_table["x", "Estimate"]
		p_value <- coef_table["x", "Pr(>|t|)"]
        t_value <- coef_table["x", "t value"]

		
		# Extract confidence intervals
		ci <- confint(fit, parm = "x", level = 0.95)
		ci_lower <- ci[1]
		ci_upper <- ci[2]
		
		# Return results
		return(data.frame(
			feature = feature,
			beta = beta,
			p.value = p_value,
			ci_lower = ci_lower,
			ci_upper = ci_upper,
            t_value = t_value,
			stringsAsFactors = FALSE
		))
        
	})
	
	# Combine results and return
	results_df <- bind_rows(results_list) %>%
        mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
        arrange(p.value)

	return(results_df)
}


test_feature_regression <- function(feature_x_cell, cell_values, ncores = 1) {
    feature_x_cell = as.matrix(feature_x_cell)
	if (is.null(names(cell_values))) stop("cell_values must be a named numeric vector")
	cells <- Reduce(intersect, list(colnames(feature_x_cell), names(cell_values)))
    message(glue::glue('{length(cells)} cells left after intersecting'))

	all_features <- rownames(feature_x_cell)

	results_list <- mclapply(all_features, mc.cores = ncores, function(feature) {

        df <- data.frame(
            y = feature_x_cell[feature, cells, drop = TRUE],
            x = cell_values[cells]
        )
        
		fit <- lm(y ~ x, data = df)
		
		# Extract summary statistics
		summary_fit <- summary(fit)
		coef_table <- summary_fit$coefficients
		
		# Extract beta coefficient and p-value for the 'x' variable
		beta <- coef_table["x", "Estimate"]
		p_value <- coef_table["x", "Pr(>|t|)"]
		
		# Extract confidence intervals
		ci <- confint(fit, parm = "x", level = 0.95)
		ci_lower <- ci[1]
		ci_upper <- ci[2]
		
		# Return results
		return(data.frame(
			feature = feature,
			beta = beta,
			p.value = p_value,
			ci_lower = ci_lower,
			ci_upper = ci_upper,
			stringsAsFactors = FALSE
		))
        
	})
	
	# Combine results and return
	results_df <- bind_rows(results_list) %>%
        mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
        arrange(p.value)

	return(results_df)
}


plot_volcano = function(test_dat, title = '', highlight = NULL) {

   if (!is.null(highlight)) {
        test_dat_label = test_dat %>% filter(p.adj < 0.05 | feature %in% highlight)
    } else {
        test_dat_label = test_dat %>% filter(p.adj < 0.05)
    }

    ggplot(
        test_dat,
        aes(x = beta, y = -log10(p.adj))
    ) +
    geom_vline(xintercept = 0, color = 'blue', linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), color = 'blue', linetype = 'dashed') +
    geom_point(
        aes(color = p.adj < 0.05)
    ) +
    geom_text_repel(
        data = test_dat_label,
        aes(label = feature),
        size = 2,
        max.overlaps = 100
    ) +
    scale_color_manual(values = c('gray', 'firebrick')) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    ggtitle(title)
}


plot_forest_regression <- function(results_df, title = "Forest Plot: Feature Associations") {
    # results_df: output from test_feature_regression_batch function
    # top_n: number of top features to plot (by p-value)

    # Filter and prepare data
    plot_data <- results_df %>%
        mutate(
            feature = factor(feature, levels = rev(feature)),
            significant = p.adj < 0.05,
            color_group = case_when(
                p.adj < 0.001 ~ "Q < 0.001",
                p.adj < 0.01 ~ "Q < 0.01",
                p.adj < 0.05 ~ "Q < 0.05",
                TRUE ~ "NS"
            )
        )

    # Create forest plot
    p <- ggplot(plot_data, aes(x = beta, y = feature)) +
        geom_vline(xintercept = 0, linetype = "solid", color = "gray50", alpha = 0.7) +
        geom_point(aes(color = color_group), size = 2) +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = color_group),
                      width = 0.2) +
        scale_color_manual(
            values = c("Q < 0.001" = "#d73027",
                      "Q < 0.01" = "#f46d43",
                      "Q < 0.05" = "#fdae61",
                      "NS" = "#abd9e9"),
            name = ""
        ) +
        labs(
            title = title,
            x = "Beta Coefficient (95% CI)",
            y = "Feature"
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(color = 'gray', linewidth = 0.25, linetype = 'dashed')
        )

    return(p)
}



matrix_ccc <- function(G, E, diag = TRUE) {

	common_rows <- intersect(rownames(G), rownames(E))
	common_cols <- intersect(colnames(G), colnames(E))

	G <- G[common_rows, common_cols, drop=FALSE]
	E <- E[common_rows, common_cols, drop=FALSE]

	ix <- upper.tri(G, diag = diag)
	g  <- G[ix]; e <- E[ix]
	mg <- mean(g); me <- mean(e)
	vg <- stats::var(g); ve <- stats::var(e)
	cge <- stats::cov(g, e)
	(2 * cge) / (vg + ve + (mg - me)^2)
}


matrix_cor_test <- function(A, B, nperm = 0, diag = TRUE) {
	# 1) Basic checks
	if (!is.matrix(A) || !is.matrix(B)) {
		stop("Both A and B must be matrices")
	}
	
	if (is.null(rownames(A)) || is.null(rownames(B)) ||
	    is.null(colnames(A)) || is.null(colnames(B))) {
		stop("Both matrices must have rownames and colnames")
	}

	# 2) Find common rows/columns (preserve A's order)
	common_rows <- intersect(rownames(A), rownames(B))
	common_cols <- intersect(colnames(A), colnames(B))

	if (length(common_rows)==0 || length(common_cols)==0) {
		stop("No overlapping rownames or colnames between A and B")
	}

	# 3) Subset & align
	A2 <- A[common_rows, common_cols, drop=FALSE]
	B2 <- B[common_rows, common_cols, drop=FALSE]

	# 4) Compute correlation on upper triangle and diagonal only
	upper_tri_idx <- upper.tri(A2, diag = diag)
	test <- cor.test(
		as.vector(A2[upper_tri_idx]),
		as.vector(B2[upper_tri_idx])
	)

    # test = vegan::mantel(A2, B2, permutations = nperm, method = "pearson")

	return(test)
}

recall_at_prec <- function(D, target_prec) {
    df <- D[order(D$prec), c("conf","prec","recall")]
    df <- df[is.finite(df$prec) & is.finite(df$recall), ]
    df <- df[df$prec > 0 & df$recall > 0, ]
    stopifnot(nrow(df) >= 2)

    i_hi <- which(df$prec >= target_prec)[1]
    if (is.na(i_hi)) return(list(mode="above_max", row=df[nrow(df), ]))
    if (i_hi == 1)   return(list(mode="below_min", row=df[1, ]))

    lo <- df[i_hi - 1, ]
    hi <- df[i_hi, ]
    t <- (target_prec - lo$prec) / (hi$prec - lo$prec)

    data.frame(
      recall_hat = lo$recall + t * (hi$recall - lo$recall),
      conf_hat   = lo$conf   + t * (hi$conf   - lo$conf)
    ) %>%
    mutate(
        recall_hat = signif(recall_hat, 3)
    )
  }