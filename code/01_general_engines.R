# 01_general_engines.R
# Purpose: High-Dimensional Analysis Functions (Agnostic / General Purpose)
# Architecture: Topic-based modules (Wrappers followed by associated Workers).

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(patchwork)
library(glmnet)
library(umap)
library(cluster)
library(dbscan)
library(igraph)
library(missForest)
library(circlize)
library(pheatmap)
library(ComplexHeatmap)

# WRAPPERS 1: Dimensionality Reduction & Utilities ----------------------------------

# 1.1 PCA Pipeline Wrapper
# Description:
#   Executes a complete Principal Component Analysis pipeline. It filters data,
#   calculates PCA, and saves a suite of dynamically scaled visualizations to PNG,
#   returning the ggplot objects and their direct file paths.
# Args:
#   omics_mat: Numeric matrix (features in rows, subjects in columns).
#   clin_df: Aligned clinical metadata dataframe.
#   pca_conf: List of settings (title, subset_to_contrast, color_vars).
#   pca_id: String identifier for saving outputs.
#   base_output_dir: Path for saving results.
#   project_colors_func: Function to map cohort specific hex codes.
# Returns:
#   Nested list with 'plots' (ggplot objects) and 'paths' (string file paths).
wrap_pca_pipeline <- function(omics_mat, clin_df, pca_conf, pca_id, base_output_dir, project_colors_func) {

  # Setup Directories
  pca_base_dir <- file.path(base_output_dir, "PCA")
  pca_sub_dir <- file.path(pca_base_dir, pca_id)
  if (!dir.exists(pca_sub_dir)) dir.create(pca_sub_dir, recursive = TRUE)

  # Data Filtering & Subset
  if (!is.null(pca_conf$subset_to_contrast)) {
    target_col <- pca_conf$subset_to_contrast$group_col
    target_lvls <- pca_conf$subset_to_contrast$target_groups
    clin_df <- clin_df %>% filter(!!sym(target_col) %in% target_lvls)
  }

  valid_barcodes <- intersect(colnames(omics_mat), clin_df$Subject_ID)
  clin_df <- clin_df %>% filter(Subject_ID %in% valid_barcodes)
  mat_sub <- omics_mat[, clin_df$Subject_ID, drop = FALSE]

  # PCA Calculation
  pca_res <- FactoMineR::PCA(t(mat_sub), scale.unit = TRUE, graph = FALSE)
  plots <- list()
  paths <- list() # Store paths for HTML rendering

  # --- AESTHETIC REFACTOR: TITLE WRAPPING ---
  # We wrap the main title at ~55 characters so it never overflows.
  wrapped_main_title <- stringr::str_wrap(pca_conf$title, width = 55)

  # Generate Individual Plots
  plots[["scree"]] <- fviz_eig(pca_res, addlabels = TRUE) +
    theme_project_base() +
    labs(title = sprintf("%s: Variance Explained (Scree)", wrapped_main_title))

  plots[["loadings"]] <- fviz_pca_var(pca_res, select.var = list(contrib = 30),
                                      col.var = "contrib",
                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                      repel = TRUE) +
    theme_project_base() +
    guides(colour = guide_colorbar(barwidth = unit(8, "cm"), barheight = unit(0.3, "cm"))) +
    labs(title = sprintf("%s: PCA Loadings (Top 30)", wrapped_main_title))

  plots[["pc1_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 1, top = 10, fill = "grey30", color = "grey30") +
    theme_project_base() + coord_flip() + labs(title = "Top PC1 Drivers")

  plots[["pc2_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 2, top = 10, fill = "grey30", color = "grey30") +
    theme_project_base() + coord_flip() + labs(title = "Top PC2 Drivers")

  # Variable Overlays (Individual and Merged)
  overlay_plots_for_grid <- list()
  for (c_var in pca_conf$color_vars) {
    if (!(c_var %in% colnames(clin_df))) next

    # --- AESTHETIC REFACTOR: SPLIT TITLE/SUBTITLE ---
    # Failsafe: Replicate design decision for interactive viewing
    interact_main_title <- stringr::str_wrap(pca_conf$title, width = 50)
    interact_labs_list <- labs(title = interact_main_title, subtitle = paste("Colored by:", c_var))

    clean_vec <- clin_df[[c_var]][!is.na(clin_df[[c_var]])]
    if (length(clean_vec) == 0) next

    # A. CONTINUOUS VARIABLE LOGIC (Reverted to Viridis from config)
    if (is.numeric(clin_df[[c_var]])) {
      cont_pal <- if(exists("FALLBACK_CONT_PALETTE")) FALLBACK_CONT_PALETTE else "viridis"

      p <- fviz_pca_ind(pca_res, label = "none", col.ind = clin_df[[c_var]], legend.title = c_var) +
        scale_color_viridis_c(option = cont_pal, na.value = "grey85") +
        theme_project_base() + interact_labs_list

      # B. CATEGORICAL VARIABLE LOGIC
    } else {
      clin_df[[c_var]] <- as.factor(clin_df[[c_var]])
      lvl_names <- levels(clin_df[[c_var]])

      mapped_colors <- project_colors_func(lvl_names)

      if (all(mapped_colors == "#333333")) {
        safe_pal <- if(exists("FALLBACK_CAT_PALETTE")) FALLBACK_CAT_PALETTE else c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#999999")
        if (length(lvl_names) <= length(safe_pal)) {
          mapped_colors <- safe_pal[1:length(lvl_names)]
        } else {
          mapped_colors <- colorRampPalette(safe_pal)(length(lvl_names))
        }
      }

      p <- fviz_pca_ind(pca_res, label = "none", habillage = clin_df[[c_var]],
                        palette = mapped_colors,
                        addEllipses = TRUE) +
        theme_project_base() + interact_labs_list
    }

    plots[[paste0("PCA_Plot_", c_var)]] <- p

    # Grid Title Setup (Cleaner than interactive viewer)
    overlay_plots_for_grid[[c_var]] <- p + labs(title = wrapped_main_title, subtitle = paste("Colored by:", c_var))
  }

  # Automated Saving (With filename sanitization)
  for (p_name in names(plots)) {
    # Sanitize to prevent '%' crashes
    safe_p_name <- gsub("%", "pct", p_name)
    safe_p_name <- gsub("[^A-Za-z0-9_.-]", "_", safe_p_name)
    safe_p_name <- gsub("_+", "_", safe_p_name)

    file_path <- file.path(pca_sub_dir, paste0(pca_id, "_", safe_p_name, ".png"))
    ggsave(filename = file_path, plot = plots[[p_name]], width = 7, height = 6, dpi = 300, bg = "white")
    paths[[p_name]] <- file_path # Store path
  }

  # Stitch & Save Merged Drivers
  driver_grid <- (plots[["pc1_drivers"]] | plots[["pc2_drivers"]]) +
    plot_annotation(title = sprintf("%s: PCA Variance Drivers", wrapped_main_title),
                    tag_levels = 'A', theme = theme(plot.title = element_text(size = 16, face = "bold")))

  driver_path <- file.path(pca_sub_dir, paste0(pca_id, "_Merged_Drivers.png"))
  ggsave(driver_path, driver_grid, width = 12, height = 6, dpi = 300, bg="white")

  plots[["merged_drivers"]] <- driver_grid
  paths[["merged_drivers"]] <- driver_path

  # Stitch & Save Merged Overlays (Dynamic height!)
  if (length(overlay_plots_for_grid) > 0) {
    ncol_val <- if(length(overlay_plots_for_grid) > 2) 2 else 1
    master_scatter_grid <- wrap_plots(overlay_plots_for_grid, ncol = ncol_val) +
      plot_annotation(title = wrapped_main_title, subtitle = "Hierarchical Grid: Metadata Colored Overlays",
                      tag_levels = 'A', theme = theme(plot.title = element_text(size = 18, face = "bold")))

    overlay_path <- file.path(pca_sub_dir, paste0(pca_id, "_Merged_Scatter_Overlays.png"))
    # Calculation handles odd number of plots gracefully
    ggsave(overlay_path, master_scatter_grid,
           width = 14, height = (6 * ceiling(length(overlay_plots_for_grid)/2)), dpi = 300, bg="white")

    plots[["merged_overlays"]] <- master_scatter_grid
    paths[["merged_overlays"]] <- overlay_path
  }

  return(list(plots = plots, paths = paths))
}

# 1.2 FAMD PIPELINE (Factor Analysis of Mixed Data)
wrap_famd_pipeline <- function(df, id_col, famd_conf, famd_id, base_output_dir, project_colors_func) {

  famd_sub_dir <- file.path(base_output_dir, "FAMD", famd_id)
  if (!dir.exists(famd_sub_dir)) dir.create(famd_sub_dir, recursive = TRUE)

  clin_data <- as.data.frame(df %>% select(-all_of(id_col)) %>% mutate_if(is.character, as.factor))
  rownames(clin_data) <- df[[id_col]]

  # FACTOMINER SAFETY PATCH: Ensure globally unique factor levels
  for (col in colnames(clin_data)) {
    if (is.factor(clin_data[[col]])) levels(clin_data[[col]]) <- paste(col, levels(clin_data[[col]]), sep = "_")
  }

  max_possible_dims <- min(nrow(clin_data) - 1, ncol(clin_data) * 2)
  famd_res <- FAMD(clin_data, graph = FALSE, ncp = min(max_possible_dims, famd_conf$max_dims))

  plots <- list(); paths <- list()
  wrapped_title <- stringr::str_wrap(famd_conf$title, width = 55)

  t_dim <- if(!is.null(famd_conf$top_dim)) min(famd_conf$top_dim, ncol(famd_res$ind$coord)) else 5
  cum_var <- famd_res$eig[t_dim, 3]

  plots[["scree"]] <- fviz_screeplot(famd_res, addlabels = TRUE) +
    theme_project_base() +
    geom_vline(xintercept = t_dim + 0.5, linetype = "dashed", color = "firebrick", linewidth = 1) +
    annotate("text", x = t_dim + 0.3, y = max(famd_res$eig[, 2]) * 0.8, label = sprintf("Cutoff (Dim %d)\nIncluded Var: %.1f%%", t_dim, cum_var), color = "firebrick", fontface = "bold", hjust = 0) +
    labs(title = sprintf("%s: Variance Explained", wrapped_title))

  plots[["var_contrib"]] <- fviz_famd_var(famd_res, "var", axes = c(1, 2), repel = TRUE, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
    theme_project_base() + labs(title = sprintf("%s: Variable Contributions (Dim 1-2)", wrapped_title))

  # Top Drivers Heatmap
  contrib_mat <- famd_res$var$contrib[, 1:t_dim, drop = FALSE]
  top_vars <- unique(as.vector(sapply(1:t_dim, function(i) head(rownames(contrib_mat)[order(contrib_mat[, i], decreasing = TRUE)], 5))))
  plot_mat <- contrib_mat[top_vars, , drop = FALSE]
  col_fun <- colorRamp2(c(0, max(plot_mat)), c("white", "firebrick4"))

  ht <- Heatmap(plot_mat, name = "Contribution (%)", column_title = sprintf("Top Drivers across First %d FAMD Dimensions", t_dim), col = col_fun, cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = TRUE, row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"), rect_gp = grid::gpar(col = "grey85", lwd = 0.5))
  plots[["drivers_heatmap"]] <- ht

  # Individual Overlays (Pure ggplot to bypass factoextra bug)
  coord_df <- as.data.frame(famd_res$ind$coord[, 1:2])
  colnames(coord_df) <- c("Dim1", "Dim2")
  coord_df[[id_col]] <- rownames(coord_df)
  plot_df <- coord_df %>% left_join(df %>% select(all_of(id_col), any_of(famd_conf$color_vars)), by = id_col)

  for (c_var in famd_conf$color_vars) {
    if (!(c_var %in% colnames(plot_df))) next
    clean_vec <- plot_df[[c_var]]
    if(all(is.na(clean_vec))) next

    if (is.numeric(clean_vec)) {
      p <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = .data[[c_var]])) + geom_point(size = 3, alpha = 0.8) + scale_color_viridis_c(na.value = "grey85") + theme_project_base() + labs(title = wrapped_title, subtitle = paste("Colored by:", c_var), x = "Dimension 1", y = "Dimension 2")
    } else {
      lvl_names <- levels(as.factor(clean_vec)); mapped_colors <- project_colors_func(lvl_names)
      if (all(mapped_colors == "#333333")) { safe_pal <- if(exists("FALLBACK_CAT_PALETTE")) FALLBACK_CAT_PALETTE else c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2"); mapped_colors <- colorRampPalette(safe_pal)(length(lvl_names)); names(mapped_colors) <- lvl_names }

      p <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = as.factor(.data[[c_var]]))) + geom_point(size = 3, alpha = 0.8) + stat_ellipse(type = "norm", linetype = 2, show.legend = FALSE, alpha = 0.5) + scale_color_manual(values = mapped_colors) + theme_project_base() + labs(title = wrapped_title, subtitle = paste("Colored by:", c_var), x = "Dimension 1", y = "Dimension 2", color = c_var)
    }

    # Save individually instead of stitching
    safe_c_var <- gsub("[^A-Za-z0-9_.-]", "_", c_var)
    ggsave(file.path(famd_sub_dir, paste0(famd_id, "_Overlay_", safe_c_var, ".png")), p, width = 7, height = 6, bg="white")

    # Store individually in the plots list
    plots[[paste0("overlay_", c_var)]] <- p
  }

  ggsave(file.path(famd_sub_dir, paste0(famd_id, "_Scree.png")), plots[["scree"]], width = 8, height = 6, bg="white")
  ggsave(file.path(famd_sub_dir, paste0(famd_id, "_VarContrib.png")), plots[["var_contrib"]], width = 10, height = 8, bg="white")
  png(file.path(famd_sub_dir, paste0(famd_id, "_Top_Drivers_Heatmap.png")), width=8, height=8, units="in", res=300); draw(ht); dev.off()

  # Return Subsetted Matrix
  famd_matrix <- t(as.data.frame(famd_res$ind$coord[, 1:t_dim, drop = FALSE]))
  return(list(status = "success", matrix = famd_matrix, plots = plots, famd_obj = famd_res))
}

# 1.3 Matrix Factory Engine
# Description: Prepares custom matrices by applying dynamic cohort subsetting,
# feature selection (strict DEA biological thresholds or Pi-Score fallbacks),
# and confounder transformations (Pruning/Residualizing). Generates an audit table.
wrap_matrix_factory <- function(omics_mat, clin_df, conf, prep_id) {
  log_msg <- c()

  # 1. Dynamic Cohort Subsetting
  if (!is.null(conf$subset_var)) {
    var_name <- conf$subset_var
    if (!(var_name %in% colnames(clin_df))) {
      return(list(status="failed", error_msg=sprintf("Subset variable '%s' not found.", var_name)))
    }

    if (!is.null(conf$subset_groups) && length(conf$subset_groups) > 0) {
      clin_sub <- clin_df %>% filter(!!sym(var_name) %in% conf$subset_groups) %>% droplevels()
      log_msg <- c(log_msg, sprintf("  -> Cohort Subset: Filtered '%s' to [%s].", var_name, paste(conf$subset_groups, collapse=", ")))
    } else if (!is.null(conf$subset_inclusion_condition)) {
      op <- conf$subset_inclusion_condition$operator
      val <- as.numeric(conf$subset_inclusion_condition$value)

      clin_sub <- switch(op,
                         ">"  = clin_df %>% filter(!!sym(var_name) > val),
                         "<"  = clin_df %>% filter(!!sym(var_name) < val),
                         ">=" = clin_df %>% filter(!!sym(var_name) >= val),
                         "<=" = clin_df %>% filter(!!sym(var_name) <= val),
                         "==" = clin_df %>% filter(!!sym(var_name) == val),
                         "="  = clin_df %>% filter(!!sym(var_name) == val),
                         clin_df)
      log_msg <- c(log_msg, sprintf("  -> Cohort Subset: Filtered where '%s' %s %s.", var_name, op, val))
    } else {
      clin_sub <- clin_df
    }
  } else {
    clin_sub <- clin_df
  }

  valid_ids <- intersect(colnames(omics_mat), clin_sub$Subject_ID)
  clin_sub <- clin_sub %>% filter(Subject_ID %in% valid_ids)
  mat_sub <- omics_mat[, valid_ids, drop=FALSE]
  log_msg <- c(log_msg, sprintf("  -> Cohort Survival: %d patients remaining.", length(valid_ids)))

  # --- INITIALIZE AUDIT TABLE ---
  audit_tbl <- data.frame(Feature = character(), logFC = numeric(), FDR = numeric(),
                          Transform_Metric = character(), Final_Status = character(), stringsAsFactors = FALSE)

  # 2. Base Feature Selection
  if (!is.null(conf$base_features)) {
    if (conf$base_features$method == "dea") {
      dea_path <- file.path(conf$base_features$dea_path, conf$base_features$source_file)
      if (!file.exists(dea_path)) return(list(status="failed", error_msg=sprintf("DEA file not found: %s", dea_path)))

      dea_res <- readRDS(dea_path) %>%
        mutate(pi_score = abs(logFC) * -log10(adj.P.Val)) %>%
        arrange(desc(pi_score))

      top_n <- conf$base_features$top_n
      fdr_cut <- conf$base_features$fdr_cutoff
      fc_cut <- conf$base_features$fc_cutoff
      min_feats <- if(!is.null(conf$base_features$min_features)) conf$base_features$min_features else 30

      if (!is.null(fdr_cut) && !is.null(fc_cut)) {
        passing_df <- dea_res %>% filter(adj.P.Val < fdr_cut & abs(logFC) > fc_cut)
        n_pass <- nrow(passing_df)

        if (n_pass < min_feats) {
          log_msg <- c(log_msg, sprintf("  -> WARNING: Only %d proteins passed thresholds. Falling back to top %d by Pi-Score.", n_pass, min_feats))
          top_feats <- dea_res %>% head(min_feats) %>% pull(Feature)
        } else if (n_pass > top_n) {
          log_msg <- c(log_msg, sprintf("  -> WARNING: %d proteins passed. Capping at top %d strictly from passing group.", n_pass, top_n))
          top_feats <- passing_df %>% head(top_n) %>% pull(Feature)
        } else {
          log_msg <- c(log_msg, sprintf("  -> Base Signature: %d proteins passed strict biological thresholds.", n_pass))
          top_feats <- passing_df %>% pull(Feature)
        }
      } else {
        log_msg <- c(log_msg, sprintf("  -> Base Signature: Top %d proteins extracted purely by Pi-Score.", top_n))
        top_feats <- dea_res %>% head(top_n) %>% pull(Feature)
      }

      valid_features <- intersect(rownames(mat_sub), top_feats)
      mat_sub <- mat_sub[valid_features, , drop=FALSE]

      # Populate Audit Table with the candidates we initially selected
      audit_tbl <- dea_res %>% filter(Feature %in% valid_features) %>%
        select(Feature, logFC, FDR = adj.P.Val) %>%
        mutate(Transform_Metric = "-", Final_Status = "Included")

    } else if (conf$base_features$method == "variance") {
      feat_vars <- apply(mat_sub, 1, var, na.rm = TRUE)
      top_feats <- names(sort(feat_vars, decreasing = TRUE))[1:min(conf$base_features$top_n, length(feat_vars))]
      valid_features <- intersect(rownames(mat_sub), top_feats)
      mat_sub <- mat_sub[valid_features, , drop=FALSE]

      # Dummy audit table since DEA isn't used
      audit_tbl <- data.frame(Feature = valid_features, logFC = NA, FDR = NA,
                              Transform_Metric = "-", Final_Status = "Included", stringsAsFactors = FALSE)
      log_msg <- c(log_msg, sprintf("  -> Base Signature: %d proteins selected via Local Variance.", length(valid_features)))
    }
  }

  # 3. Transformations (Age Confounder Correction)
  if (!is.null(conf$transform)) {

    # --- PRUNING LOGIC ---
    if (conf$transform$method == "prune") {
      cov_var <- conf$transform$covariate
      p_thresh <- if(!is.null(conf$transform$p_threshold)) conf$transform$p_threshold else 0.05

      if (cov_var %in% colnames(clin_sub)) {
        cov_vec <- as.numeric(clin_sub[[cov_var]])
        drop_feats <- c()

        for (f in rownames(mat_sub)) {
          # Use cor.test to get both rho and p-value
          ct <- tryCatch(cor.test(mat_sub[f, ], cov_vec, method="spearman", exact=FALSE), error=function(e) NULL)
          if (!is.null(ct)) {
            rho <- unname(ct$estimate)
            p_val <- ct$p.value

            # Record in the audit table
            audit_tbl$Transform_Metric[audit_tbl$Feature == f] <- sprintf("rho=%.2f, p=%.3f", rho, p_val)

            # STRICT CHECK: Drop only if strongly correlated AND statistically significant
            if (abs(rho) > conf$transform$threshold && p_val < p_thresh) {
              drop_feats <- c(drop_feats, f)
              audit_tbl$Final_Status[audit_tbl$Feature == f] <- "Excluded (Pruned)"
            }
          }
        }
        keep_feats <- setdiff(rownames(mat_sub), drop_feats)
        mat_sub <- mat_sub[keep_feats, , drop=FALSE]
        log_msg <- c(log_msg, sprintf("  -> Transformation (Pruned): %d proteins removed (|rho| > %.2f & p < %s). **%d proteins remain.**", length(drop_feats), conf$transform$threshold, p_thresh, nrow(mat_sub)))
      }

      # --- RESIDUALIZATION LOGIC ---
    } else if (conf$transform$method == "residualize") {
      cov_vars <- conf$transform$covariates
      existing_covs <- intersect(cov_vars, colnames(clin_sub))
      if (length(existing_covs) > 0) {
        cov_df <- clin_sub %>% select(all_of(existing_covs)) %>% mutate(across(everything(), as.numeric))
        complete_cases <- complete.cases(cov_df)

        if (sum(!complete_cases) > 0) {
          cov_df <- cov_df[complete_cases, , drop=FALSE]
          mat_sub <- mat_sub[, complete_cases, drop=FALSE]
          log_msg <- c(log_msg, sprintf("  -> Warning: Dropped %d patients due to missing covariate data.", sum(!complete_cases)))
        }

        # Calculate exact linear coefficients before residualizing
        design_mat <- model.matrix(~ ., data = cov_df)
        fit <- limma::lmFit(mat_sub, design_mat)

        # Format the regression coefficients for the audit table
        for (f in rownames(mat_sub)) {
          coef_str <- paste(sapply(existing_covs, function(cv) {
            beta <- fit$coefficients[f, cv]
            sprintf("%s(beta=%.3f)", cv, beta)
          }), collapse=" | ")

          audit_tbl$Transform_Metric[audit_tbl$Feature == f] <- coef_str
        }

        suppressMessages({
          mat_sub <- limma::removeBatchEffect(mat_sub, covariates = as.matrix(cov_df))
        })
        log_msg <- c(log_msg, sprintf("  -> Transformation (Residualized): Matrix updated to residuals using %s.", paste(existing_covs, collapse=", ")))
      }
    }
  } else {
    log_msg <- c(log_msg, "  -> Transformation: None (Raw Matrix preserved).")
  }

  final_log <- paste(log_msg, collapse="\n")
  return(list(status = "success", matrix = mat_sub, log_text = final_log, audit_tbl = audit_tbl))
}

# 1.4 Global Omics Variance Filter Engine
# Description:
#   Calculates global variance for all omics features, generates a density
#   histogram, and extracts the top N highly variable features for downstream
#   unsupervised discovery.
wrap_global_variance_filter <- function(omics_mat, conf, base_output_dir) {

  # Fallback to "Feature" if omics_type isn't specified
  omics_type <- if(!is.null(conf$omics_type)) conf$omics_type else "Feature"

  # 1. Calculate variance for every feature
  feature_variances <- apply(omics_mat, 1, var, na.rm = TRUE)
  var_df <- data.frame(
    Feature = names(feature_variances),
    Variance = as.numeric(feature_variances)
  ) %>% arrange(desc(Variance))

  # 2. Determine threshold safely
  top_n <- min(conf$top_n_features, nrow(var_df))
  cutoff_val <- var_df$Variance[top_n]
  high_var_feats <- var_df$Feature[1:top_n]

  # 3. Create Plot
  p_var <- ggplot(var_df, aes(x = Variance)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey70", color = "white", alpha = 0.8) +
    geom_density(color = "#2980b9", linewidth = 1) +
    geom_vline(xintercept = cutoff_val, color = "#c0392b", linetype = "dashed", linewidth = 1) +
    annotate("text", x = cutoff_val, y = Inf,
             label = sprintf(" Cutoff: Top %d %ss\n (Var > %.2f)", top_n, omics_type, cutoff_val),
             hjust = -0.1, vjust = 1.5, color = "#c0392b", fontface = "bold") +
    theme_project_base() +
    labs(
      title = conf$title,
      subtitle = sprintf("Total %ss: %d | Selected: %d", omics_type, nrow(var_df), length(high_var_feats)),
      x = sprintf("%s Variance", omics_type),
      y = "Density"
    )

  # 4. Save Data & Plot
  safe_title <- gsub("[^A-Za-z0-9_.-]", "_", conf$title)
  plot_path <- file.path(base_output_dir, paste0("Global_Variance_", safe_title, ".png"))
  rds_path <- file.path(base_output_dir, "high_var_features.rds")

  ggsave(plot_path, p_var, width = 8, height = 6, dpi = 300, bg = "white")
  saveRDS(high_var_feats, rds_path)

  return(list(status = "success", plot = p_var, features = high_var_feats))
}

# 1.5 MISSING DATA IMPUTATION: RANDOM FOREST (missForest)
wrap_rf_imputation <- function(df, id_col = "Subject_ID", maxiter = 10, ntree = 100) {
  library(missForest)

  # 1. Isolate clinical data from the ID column
  clin_data <- df %>% select(-all_of(id_col))

  # 2. Audit: Count missingness before imputation
  missing_counts <- colSums(is.na(clin_data))
  audit_df <- data.frame(
    Variable = names(missing_counts),
    Missing_N = missing_counts,
    Missing_Pct = round((missing_counts / nrow(clin_data)) * 100, 1)
  ) %>% filter(Missing_N > 0) %>% arrange(desc(Missing_N))

  if (nrow(audit_df) == 0) {
    return(list(status = "success", imputed_df = df, audit = audit_df, error_metrics = NULL))
  }

  # 3. Type Enforcement (missForest strictly requires factors for categorical)
  clin_data <- clin_data %>% mutate_if(is.character, as.factor)

  # 4. Imputation Engine
  # Set seed inside for reproducibility of the random forest
  set.seed(42)
  rf_res <- missForest(as.data.frame(clin_data), maxiter = maxiter, ntree = ntree, verbose = FALSE)

  # 5. Reconstruct DataFrame
  imputed_df <- bind_cols(df %>% select(all_of(id_col)), as_tibble(rf_res$ximp))

  return(list(
    status = "success",
    imputed_df = imputed_df,
    audit = audit_df,
    error_metrics = rf_res$OOBerror # Normalized Root Mean Squared Error (NRMSE) & Proportion of Falsely Classified (PFC)
  ))
}

# 1.6 GOWER DISTANCE FACTORY (For Mixed Data)
wrap_gower_distance <- function(df, id_col = "Subject_ID") {
  library(cluster)

  clin_data <- df %>% select(-all_of(id_col)) %>% mutate_if(is.character, as.factor)

  # Calculate Gower Distance
  gower_dist <- cluster::daisy(clin_data, metric = "gower")

  # Convert to standard matrix for pipeline compatibility, reattaching rownames
  gower_mat <- as.matrix(gower_dist)
  rownames(gower_mat) <- df[[id_col]]
  colnames(gower_mat) <- df[[id_col]]

  return(list(
    dist_obj = gower_dist,
    matrix = gower_mat
  ))
}



# WRAPPERS 2: Clustering & Predictive Modeling -----------------------------------

# 2.1 K-Optimization (Pre-Flight) Pipeline Wrapper
# Description:
#   Algorithm-Aware Optimization. Sweeps 'k' for classic algorithms (Silhouette/WSS),
#   sweeps 'minPts' for HDBSCAN (Noise vs Clusters), and sweeps 'resolution' for
#   Louvain (Modularity).
wrap_k_optimization_pipeline <- function(omics_mat, conf, opt_id, base_output_dir) {
  library(ggplot2)
  library(patchwork)
  library(cluster)

  results <- list(status = "success", plots = list(), paths = list())

  opt_dir <- file.path(base_output_dir, "K_Optimization", opt_id)
  if (!dir.exists(opt_dir)) dir.create(opt_dir, recursive = TRUE)

  # Check if we are handling a pre-computed distance matrix (e.g., Gower)
  is_dist <- if (!is.null(conf$is_dist_matrix)) conf$is_dist_matrix else FALSE
  patient_mat <- if (isTRUE(is_dist)) omics_mat else t(omics_mat)

  for (algo_name in names(conf$algorithms)) {
    algo_conf <- conf$algorithms[[algo_name]]
    if (!isTRUE(algo_conf$run)) next

    combined_plot <- NULL

    # --- A. HDBSCAN Optimization Sweep ---
    if (algo_name == "hdbscan" && !isTRUE(is_dist)) {
      library(dbscan)
      minPts_vals <- if(!is.null(algo_conf$minPts_range)) algo_conf$minPts_range else 3:15
      sweep_df <- data.frame(minPts = minPts_vals, Clusters = 0, Noise_Pct = 0)

      for (i in seq_along(minPts_vals)) {
        hd <- dbscan::hdbscan(patient_mat, minPts = minPts_vals[i])
        sweep_df$Clusters[i] <- max(hd$cluster)
        sweep_df$Noise_Pct[i] <- sum(hd$cluster == 0) / length(hd$cluster)
      }

      max_c <- max(sweep_df$Clusters, 1)
      combined_plot <- ggplot(sweep_df, aes(x = minPts)) +
        geom_line(aes(y = Clusters, color = "Clusters Found"), linewidth = 1) +
        geom_point(aes(y = Clusters, color = "Clusters Found"), size = 3) +
        geom_line(aes(y = Noise_Pct * max_c, color = "% Noise"), linewidth = 1, linetype="dashed") +
        scale_y_continuous(name = "Number of Clusters", sec.axis = sec_axis(~ . / max_c, name = "% Noise", labels = scales::percent)) +
        scale_color_manual(values = c("Clusters Found" = "#2980b9", "% Noise" = "#c0392b")) +
        theme_project_base() + theme(legend.position = "bottom", legend.title = element_blank()) +
        labs(title = "HDBSCAN Optimization", subtitle = "Sweep: minPts (Look for cluster stability with low noise)", x = "minPts")

      # --- B. LOUVAIN Optimization Sweep ---
    } else if (algo_name == "louvain" && !isTRUE(is_dist)) {
      library(dbscan)
      library(igraph)
      res_vals <- if(!is.null(algo_conf$resolution_range)) algo_conf$resolution_range else seq(0.1, 2.0, by=0.1)
      knn_k <- if(!is.null(algo_conf$k_neighbors)) algo_conf$k_neighbors else 10

      knn_obj <- dbscan::kNN(patient_mat, k = knn_k)
      edges <- data.frame(from = rep(1:nrow(patient_mat), each = knn_k), to = as.vector(t(knn_obj$id)), weight = 1 / (1 + as.vector(t(knn_obj$dist))))
      g <- igraph::simplify(igraph::graph_from_data_frame(edges, directed = FALSE), remove.multiple = TRUE, remove.loops = TRUE)

      sweep_df <- data.frame(Resolution = res_vals, Modularity = 0, Clusters = 0)
      for (i in seq_along(res_vals)) {
        lc <- igraph::cluster_louvain(g, resolution = res_vals[i])
        sweep_df$Modularity[i] <- igraph::modularity(lc)
        sweep_df$Clusters[i] <- length(unique(lc$membership))
      }

      combined_plot <- ggplot(sweep_df, aes(x = Resolution)) +
        geom_line(aes(y = Modularity), color = "#27ae60", linewidth = 1) +
        geom_point(aes(y = Modularity), color = "#27ae60", size = 3) +
        geom_text(aes(y = Modularity, label = paste0(Clusters, "c")), vjust = -1, size = 3.5, color = "grey30") +
        theme_project_base() +
        labs(title = "Graph Community Optimization", subtitle = "Sweep: Resolution (Peak Modularity = Best Web Split)", x = "Resolution Parameter", y = "Modularity Score")

      # --- C. PRE-COMPUTED DISTANCE SWEEP (Gower/Custom) ---
    } else if (isTRUE(is_dist)) {
      # Skip algorithms that require continuous coordinates
      if (algo_name %in% c("hdbscan", "louvain", "kmeans")) next

      k_max <- max(conf$k_range)
      sil_scores <- numeric(k_max)
      dist_obj <- as.dist(patient_mat)

      for(k_val in 2:k_max) {
        labels <- run_clustering_worker(patient_mat, k_val, algo_name, algo_conf, is_dist_matrix = TRUE)
        num_labels <- as.numeric(as.factor(labels))
        if(length(unique(num_labels)) > 1) {
          sil_scores[k_val] <- mean(cluster::silhouette(num_labels, dist_obj)[, "sil_width"])
        }
      }

      sweep_df <- data.frame(k = 2:k_max, Silhouette = sil_scores[2:k_max])
      combined_plot <- ggplot(sweep_df, aes(x = k, y = Silhouette)) +
        geom_line(linewidth = 1, color = "#2980b9") +
        geom_point(size = 3, color = "#2980b9") +
        theme_project_base() +
        scale_x_continuous(breaks = 2:k_max) +
        labs(title = "Average Silhouette Width", subtitle = sprintf("Algorithm: %s (Distance Matrix)", toupper(algo_name)), x = "Number of Clusters (k)")

      # --- D. CLASSIC Optimization Sweep (K-Means/PAM/H-Clust on Features) ---
    } else {
      k_max <- max(conf$k_range)
      gap_boot <- if(!is.null(conf$gap_bootstraps)) conf$gap_bootstraps else 0

      FUNcluster <- function(x, k) {
        labels <- run_clustering_worker(x, k, algo_name, algo_conf, is_dist_matrix = FALSE)
        list(cluster = as.numeric(as.factor(labels)))
      }

      p_wss <- factoextra::fviz_nbclust(patient_mat, FUNcluster, method = "wss", k.max = k_max) +
        theme_project_base() + labs(title = "Elbow Method (WSS)", subtitle = sprintf("Algorithm: %s", toupper(algo_name)))

      p_sil <- factoextra::fviz_nbclust(patient_mat, FUNcluster, method = "silhouette", k.max = k_max) +
        theme_project_base() + labs(title = "Average Silhouette", subtitle = sprintf("Algorithm: %s", toupper(algo_name)))

      if (gap_boot > 0) {
        set.seed(42)
        p_gap <- factoextra::fviz_nbclust(patient_mat, FUNcluster, method = "gap_stat", k.max = k_max, nboot = gap_boot) +
          theme_project_base() + labs(title = "Gap Statistic", subtitle = sprintf("Bootstraps: %d", gap_boot))
        combined_plot <- (p_wss | p_sil | p_gap) + plot_annotation(title = sprintf("K-Opt: %s (%s)", conf$title, toupper(algo_name)))
      } else {
        combined_plot <- (p_wss | p_sil) + plot_annotation(title = sprintf("K-Opt: %s (%s)", conf$title, toupper(algo_name)))
      }
    }

    # Save & Export
    if (!is.null(combined_plot)) {
      safe_algo <- gsub("[^A-Za-z0-9_.-]", "_", algo_name)
      plot_path <- file.path(opt_dir, paste0(opt_id, "_", safe_algo, "_k_opt.png"))
      w <- if(algo_name %in% c("hdbscan", "louvain") || isTRUE(is_dist)) 8 else if(!is.null(conf$gap_bootstraps) && conf$gap_bootstraps > 0) 15 else 10
      ggsave(plot_path, combined_plot, width = w, height = 5, dpi = 300, bg = "white")

      results$plots[[algo_name]] <- combined_plot
      results$paths[[algo_name]] <- plot_path
    }
  }

  return(results)
}

# 2.2 UMAP & Multi-Algorithm Molecular Clustering Engine
# Description:
#   Performs unsupervised clustering on TRUE high-dimensional space, calculates a 2D UMAP
#   strictly for visualization, and paints the high-dimensional clusters onto the UMAP.
#   NEW: Automatically generates UMAPs colored by clinical metadata if provided.
wrap_clustering_pipeline <- function(target_matrix, clin_df, job_conf, job_id, base_output_dir, project_colors_func) {
  library(umap)
  library(ggplot2)
  library(patchwork)
  library(dplyr)

  # 1. Setup Dynamic Output Directories
  cluster_dir <- file.path(base_output_dir, "Clustering", job_id)
  if (!dir.exists(cluster_dir)) dir.create(cluster_dir, recursive = TRUE)

  # 2. Safely Extract Distance Flag (Preserves legacy script compatibility)
  is_dist <- if (!is.null(job_conf$is_dist_matrix)) job_conf$is_dist_matrix else FALSE

  # 3. Orient Matrix (Ensure Rows = Patients, Columns = Features/Dimensions)
  patient_mat <- if (isTRUE(is_dist)) target_matrix else t(target_matrix)

  # 4. Execute Algorithms & Initialize Plot List EARLY
  assignments <- list()
  plots <- list(umaps = list(), dendrograms = list()) # <--- FIX: INITIALIZED HERE

  for (algo in names(job_conf$algorithms)) {
    algo_conf <- job_conf$algorithms[[algo]]

    if (!isTRUE(algo_conf$run)) next

    # Pass the flag down to the worker
    cluster_labels <- run_clustering_worker(
      patient_matrix = patient_mat,
      k = job_conf$k,
      algo_name = algo,
      algo_conf = algo_conf,
      is_dist_matrix = is_dist
    )

    if (!is.null(cluster_labels)) {
      assignments[[algo]] <- cluster_labels

      # --- DENDROGRAM EXTRACTION ---
      tree_obj <- attr(cluster_labels, "tree")
      if (!is.null(tree_obj)) {
        if (!exists("dendrograms", where = plots)) plots$dendrograms <- list()

        # Fetch exact project colors for the branches
        lvl_names <- levels(as.factor(cluster_labels))
        mapped_colors <- project_colors_func(lvl_names)
        if (all(mapped_colors == "#333333")) mapped_colors <- NULL

        p_dend <- factoextra::fviz_dend(
          tree_obj, k = job_conf$k,
          palette = mapped_colors,
          rect = TRUE, rect_fill = TRUE,
          rect_border = mapped_colors,
          cex = 0.5, # Small text for patient barcodes
          main = sprintf("%s: Hierarchical Dendrogram (k=%s)", job_conf$title, job_conf$k)
        ) +
          theme_project_base() +
          theme(axis.text.x = element_blank()) +
          guides(size = "none") # <--- FIX 1: Destroys the black square dummy legend!

        # <--- FIX 2: Intercept the ggplot layers and forcefully inject modern linewidth!
        for (i in seq_along(p_dend$layers)) {
          if (inherits(p_dend$layers[[i]]$geom, "GeomSegment")) {
            p_dend$layers[[i]]$aes_params$linewidth <- 1 # Adjust this exact number for thickness
          }
        }

        # Export and Store
        ggsave(file.path(cluster_dir, sprintf("Dendrogram_%s.png", algo)), plot = p_dend, width = 10, height = 6, bg = "white", dpi = 300)
        plots$dendrograms[[algo]] <- p_dend
      }
    }
  }

  # Failsafe
  if (length(assignments) == 0) {
    return(list(status = "failed", error_msg = "All algorithms failed or returned NULL (likely due to algorithm/matrix incompatibility)."))
  }

  # 5. UMAP Generation (For 2D Visualization)
  n_neigh <- min(if(!is.null(job_conf$n_neighbors)) job_conf$n_neighbors else 15, nrow(patient_mat) - 1)

  if (isTRUE(is_dist)) {
    umap_res <- umap::umap(patient_mat, input = "dist", n_neighbors = n_neigh)
  } else {
    umap_res <- umap::umap(patient_mat, n_neighbors = n_neigh)
  }

  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Subject_ID <- rownames(patient_mat)

  # Bind the clinical metadata to the UMAP coordinates for coloring
  umap_df <- umap_df %>% left_join(clin_df, by = "Subject_ID")

  # 6. Generate Plots
  # (Removed the late list initialization from here)

  for (algo in names(assignments)) {

    # Isolate labels for this specific algorithm
    plot_df <- umap_df
    plot_df$Cluster <- assignments[[algo]][plot_df$Subject_ID]
    plot_df <- plot_df %>% filter(!is.na(Cluster))

    # Define Colors using the global project dictionary
    lvl_names <- levels(as.factor(plot_df$Cluster))
    mapped_colors <- project_colors_func(lvl_names)

    # A. Build the Base Cluster Plot
    p_base <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(fill = Cluster, shape = Cluster), color = "grey30", size = 3, alpha = 0.8) +
      scale_shape_manual(values = rep(21:25, length.out = length(lvl_names))) +
      scale_fill_manual(values = mapped_colors) +
      theme_project_base() +
      labs(
        title = sprintf("%s: %s", job_conf$title, toupper(algo)),
        subtitle = sprintf("k=%s | %s Space", job_conf$k, if(is_dist) "Gower Distance" else "Continuous Coordinate")
      )

    # B. Generate Clinical Overlays (If requested in config)
    if (!is.null(job_conf$color_vars)) {
      sub_plots <- list()

      for (c_var in job_conf$color_vars) {
        if (c_var %in% colnames(plot_df)) {

          if (is.numeric(plot_df[[c_var]])) {
            p_sub <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = .data[[c_var]])) +
              geom_point(size = 2, alpha = 0.8) +
              scale_color_viridis_c(na.value = "grey85") +
              theme_project_base() +
              labs(title = paste("Colored by:", c_var))
          } else {
            c_lvls <- levels(as.factor(plot_df[[c_var]]))
            c_colors <- project_colors_func(c_lvls)

            p_sub <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = as.factor(.data[[c_var]]))) +
              geom_point(size = 2, alpha = 0.8) +
              scale_color_manual(values = c_colors) +
              theme_project_base() +
              labs(title = paste("Colored by:", c_var))
          }
          sub_plots[[c_var]] <- p_sub
        }
      }

      # Stitch the base plot and overlays together using Patchwork
      if (length(sub_plots) > 0) {
        p_final <- wrap_plots(c(list(p_base), sub_plots), ncol = 2)
      } else {
        p_final <- p_base
      }

    } else {
      p_final <- p_base
    }

    # Export Dynamic PNG
    safe_algo <- gsub("[^A-Za-z0-9_.-]", "_", algo)
    ggsave(file.path(cluster_dir, sprintf("UMAP_%s.png", safe_algo)), plot = p_final, width = 12, height = 8, bg = "white")

    plots$umaps[[algo]] <- p_final
  }

  # 7. Package and Return Payload
  payload <- list(
    status = "success",
    data = list(
      assignments = assignments,
      umap_coords = umap_df,
      # Store the correct distance matrix for the downstream Diagnostics engine
      dist_matrix = if(is_dist) target_matrix else as.matrix(dist(patient_mat))
    ),
    plots = plots
  )

  return(payload)
}

# 2.3 Consensus Export & Alignment Engine
# Description: Resolves arbitrary cluster naming via Reference-Based Greedy Overlap,
# computes a multi-algorithm consensus endotype, and generates convergence diagnostics.
wrap_consensus_export <- function(clustering_payload, anchor_algo, conf, export_id, base_output_dir) {

  export_dir <- file.path(base_output_dir, "Exports")
  if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  assignments <- clustering_payload$data$assignments
  dist_mat <- clustering_payload$data$dist_matrix
  umap_df <- clustering_payload$data$umap_df

  if (is.null(assignments) || !(anchor_algo %in% names(assignments))) {
    return(list(status = "failed", error_msg = "Invalid payload or anchor algorithm missing."))
  }

  # 1. Setup Base Dataframe (Anchor is truth)
  anchor_labels <- assignments[[anchor_algo]]
  master_df <- data.frame(Subject_ID = names(anchor_labels), stringsAsFactors = FALSE)
  master_df[[paste0(toupper(anchor_algo), "_Cluster")]] <- as.character(anchor_labels)

  num_anchor <- as.numeric(as.factor(anchor_labels))
  sil_anchor <- silhouette(num_anchor, dist_mat)
  master_df[[paste0(toupper(anchor_algo), "_Sil")]] <- round(sil_anchor[, "sil_width"], 3)

  aligned_labels_list <- list()
  aligned_labels_list[[toupper(anchor_algo)]] <- as.character(anchor_labels)

  # 2. Greedy Alignment Algorithm
  for (algo in names(assignments)) {
    if (algo == anchor_algo) next

    target_labels <- as.character(assignments[[algo]])
    tbl <- table(Anchor = as.character(anchor_labels), Target = target_labels)
    mapping <- setNames(rep(NA, ncol(tbl)), colnames(tbl))

    for (i in 1:min(nrow(tbl), ncol(tbl))) {
      max_idx <- which(tbl == max(tbl), arr.ind = TRUE)[1, , drop = FALSE]
      r_name <- rownames(tbl)[max_idx[1]]
      c_name <- colnames(tbl)[max_idx[2]]

      mapping[c_name] <- r_name
      tbl[max_idx[1], ] <- -1
      tbl[, max_idx[2]] <- -1
    }

    translated_labels <- unname(mapping[target_labels])
    aligned_labels_list[[toupper(algo)]] <- translated_labels
    master_df[[paste0(toupper(algo), "_Cluster")]] <- translated_labels

    num_target <- as.numeric(as.factor(assignments[[algo]]))
    sil_target <- silhouette(num_target, dist_mat)
    master_df[[paste0(toupper(algo), "_Sil")]] <- round(sil_target[, "sil_width"], 3)
  }

  # 3. Consensus Math (Dynamically scales to N algorithms)
  algos_to_check <- names(aligned_labels_list)
  n_algos <- length(algos_to_check)

  consensus_vec <- character(nrow(master_df))
  stability_vec <- numeric(nrow(master_df))
  max_votes_vec <- numeric(nrow(master_df)) # Track raw votes for the bar chart

  for (i in 1:nrow(master_df)) {
    row_votes <- sapply(algos_to_check, function(x) master_df[[paste0(x, "_Cluster")]][i])
    vote_table <- table(row_votes)
    max_votes <- max(vote_table)
    max_votes_vec[i] <- max_votes

    if (max_votes == 1 && n_algos > 1) {
      consensus_vec[i] <- "Mixed/Transitional"
    } else {
      consensus_vec[i] <- names(vote_table)[which.max(vote_table)]
    }

    stability_vec[i] <- round(max_votes / n_algos, 3)
  }

  master_df$Consensus_Cluster <- consensus_vec
  master_df$Stability_Score <- stability_vec
  master_df$Max_Votes <- max_votes_vec

  # Apply confidence threshold flag
  if (!is.null(conf$min_stability_threshold)) {
    master_df <- master_df %>%
      mutate(Consensus_Cluster = ifelse(Stability_Score < conf$min_stability_threshold, "Mixed/Transitional", Consensus_Cluster))
  }

  # 4. Generate Diagnostics: Convergence Bar Chart
  master_df$Vote_Ratio <- paste0(master_df$Max_Votes, "/", n_algos)

  summary_df <- master_df %>%
    group_by(Vote_Ratio, Stability_Score) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Color = case_when(
        Stability_Score == 1 ~ "#27ae60", # Perfect (Green)
        Stability_Score >= conf$min_stability_threshold ~ "#f39c12", # Passed (Yellow/Orange)
        TRUE ~ "#c0392b" # Failed (Red)
      )
    ) %>% arrange(Stability_Score)

  p_bar <- ggplot(summary_df, aes(x = reorder(Vote_Ratio, Stability_Score), y = Count, fill = Color)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.85) +
    scale_fill_identity() +
    geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_project_base() +
    labs(title = "Algorithmic Convergence",
         subtitle = sprintf("Job: %s | Threshold: >= %.2f | Total Algorithms: %d", export_id, conf$min_stability_threshold, n_algos),
         x = "Agreement (Votes / Total Algorithms)", y = "Number of Patients")

  # 5. Generate Diagnostics: Consensus UMAP
  plot_df <- master_df %>% left_join(umap_df, by = "Subject_ID")

  # Dynamically map colors. Distinct clusters get colors, Mixed gets Grey.
  cluster_lvls <- setdiff(unique(plot_df$Consensus_Cluster), "Mixed/Transitional")
  safe_pal <- c("#e74c3c", "#9b59b6", "#f1c40f", "#1abc9c", "#34495e", "#e67e22")
  color_map <- setNames(safe_pal[1:length(cluster_lvls)], cluster_lvls)
  color_map["Mixed/Transitional"] <- "grey75"

  p_umap <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, fill = Consensus_Cluster)) +
    geom_point(shape = 21, color = "black", size = 3.5, alpha = 0.8) +
    scale_fill_manual(values = color_map) +
    theme_project_base() +
    labs(title = "Final Consensus Endotypes",
         subtitle = sprintf("Job: %s | Anchor Algorithm: %s | Patients: %d", export_id, toupper(anchor_algo), nrow(plot_df)),
         fill = "Consensus")

  # 6. Save and Return
  rds_path <- file.path(export_dir, paste0(export_id, ".rds"))
  csv_path <- file.path(export_dir, paste0(export_id, ".csv"))

  # Remove our temporary Max_Votes column before saving
  save_df <- master_df %>% select(-Max_Votes, -Vote_Ratio)
  saveRDS(save_df, rds_path)
  write_csv(save_df, csv_path)

  ggsave(file.path(export_dir, paste0(export_id, "_Convergence_Bar.png")), p_bar, width = 7, height = 5, bg="white")
  ggsave(file.path(export_dir, paste0(export_id, "_Consensus_UMAP.png")), p_umap, width = 8, height = 6, bg="white")

  text_log <- sprintf(
    "  -> Anchor Algorithm: %s (Automatic Alignment Applied)\n  -> Consensus Reached: %d patients classified.\n  -> Mixed/Transitional: %d patients lacked algorithmic agreement.\n  -> Export Location: %s",
    toupper(anchor_algo),
    sum(master_df$Consensus_Cluster != "Mixed/Transitional"),
    sum(master_df$Consensus_Cluster == "Mixed/Transitional"),
    export_id
  )

  return(list(status = "success", data = save_df, text_log = text_log, plots = list(convergence_bar = p_bar, consensus_umap = p_umap)))
}

# 2.4 Clinical Phenotyping Engine (ComplexHeatmap version)
wrap_clinical_phenotyping <- function(clustering_payload, clin_df, dict_df, conf, job_id, base_output_dir) {
  results <- list(status = "success", plots = list(), tables = list())
  pheno_dir <- file.path(base_output_dir, "Clinical_Phenotyping", job_id)
  if (!dir.exists(pheno_dir)) dir.create(pheno_dir, recursive = TRUE)

  for (algo_name in conf$target_algorithms) {
    labels <- clustering_payload$data$assignments[[algo_name]]
    if (is.null(labels)) next

    results$plots[[algo_name]] <- list(); results$tables[[algo_name]] <- list()
    results$plots[[algo_name]][["heatmap"]] <- list()
    results$plots[[algo_name]][["highlights"]] <- list()

    plot_df <- clin_df %>% filter(Subject_ID %in% names(labels))
    plot_df$Cluster <- as.factor(labels[plot_df$Subject_ID])
    valid_df <- plot_df %>% filter(!grepl("Noise|Mixed", Cluster)) %>% droplevels()

    # --- SUBCLUSTER FILTERING (Option A) ---
    if (!is.null(conf$subset_filter)) {
      filt_col <- conf$subset_filter$column
      filt_val <- conf$subset_filter$value
      if (filt_col %in% colnames(valid_df)) {
        valid_df <- valid_df %>% filter(.data[[filt_col]] %in% filt_val) %>% droplevels()
      }
    }

    if (length(unique(valid_df$Cluster)) < 2) next

    # --- NEW: Safely handle both flat vectors and named lists ---
    vars_list <- if (is.list(conf$variables)) conf$variables else list("Clinical_Phenotypes" = conf$variables)
    hl_list <- if (is.list(conf$highlight_variables)) conf$highlight_variables else list("Selected_Highlights" = conf$highlight_variables)

    # 1. CLINICAL TABLE (Uses all variables unlisted for a single, master summary table)
    all_vars <- unname(unlist(vars_list))
    summary_obj <- valid_df %>%
      select(any_of(all_vars), Cluster) %>%
      gtsummary::tbl_summary(by = Cluster, missing = "no") %>%
      gtsummary::add_p() %>% gtsummary::sort_p() %>% gtsummary::bold_labels()

    results$tables[[algo_name]] <- summary_obj

    # EXPORT RAW GTSUMMARY TO CSV
    raw_tbl <- gtsummary::as_tibble(summary_obj)
    readr::write_csv(raw_tbl, file.path(pheno_dir, paste0(job_id, "_", algo_name, "_Clinical_Table.csv")))

    # EXPORT RAW GTSUMMARY TO CSV
    raw_tbl <- gtsummary::as_tibble(summary_obj)
    readr::write_csv(raw_tbl, file.path(pheno_dir, paste0(job_id, "_", algo_name, "_Clinical_Table.csv")))

    # --- TARGETED HIGHLIGHT TABLES (Bypass Loop) ---
    if (!is.null(conf$highlight_variables)) {
      # FIX 1: Safely initialize the highlights list inside results$tables
      if (!("highlights" %in% names(results$tables))) results$tables$highlights <- list()
      results$tables$highlights[[algo_name]] <- list()

      for (h_name in names(conf$highlight_variables)) {
        h_vars <- conf$highlight_variables[[h_name]]
        # Filter to variables that actually exist in the data
        valid_h_vars <- intersect(h_vars, colnames(valid_df))

        if (length(valid_h_vars) > 0) {
          # FIX 2: Subset from valid_df so we use the clean "Cluster" column and exclude noise!
          hl_df <- valid_df %>% select(all_of(c("Cluster", valid_h_vars)))

          # Generate the gtsummary table showing ALL variables in this list, ignoring p_cutoff
          hl_tbl <- hl_df %>%
            gtsummary::tbl_summary(
              by = Cluster,
              missing = "ifany",
              missing_text = "Missing (NA)"
            ) %>%
            gtsummary::add_p(test = list(gtsummary::all_continuous() ~ "wilcox.test", gtsummary::all_categorical() ~ "fisher.test")) %>%
            gtsummary::bold_labels()

          # Convert to HTML and Store
          hl_html <- gtsummary::as_kable_extra(hl_tbl, format = "html", caption = sprintf("Subset Table: %s", gsub("_", " ", h_name))) %>%
            kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE, position = "left")

          # FIX 3: Save to the correct results object
          results$tables$highlights[[algo_name]][[h_name]] <- hl_html
        }
      }
    }

    # --- NEW: OUTER LOOP FOR HEATMAPS ---
    for (cat_name in names(vars_list)) {
      cat_vars <- vars_list[[cat_name]]

      # 2. ENRICHMENT DATA PREP
      enrichment_list <- list()
      mandatory_vars <- intersect(conf$must_include_vars, colnames(valid_df))
      discovery_vars <- setdiff(intersect(cat_vars, colnames(valid_df)), mandatory_vars)

      for (v in unique(c(mandatory_vars, discovery_vars))) {
        var_class <- if (v %in% dict_df$Variable) dict_df$Class[dict_df$Variable == v] else "Continuous"
        vec <- valid_df[[v]]; clust_vec <- valid_df$Cluster
        is_mandatory <- v %in% mandatory_vars

        if (grepl("Continuous|Ordinal", var_class, ignore.case = TRUE) && is.numeric(vec)) {
          p_val <- tryCatch(kruskal.test(vec ~ clust_vec)$p.value, error=function(e) NA)
          if (!is_mandatory && (is.na(p_val) || p_val >= conf$p_cutoff)) next
          g_m <- mean(vec, na.rm=T); g_sd <- sd(vec, na.rm=T)
          for (c_lvl in levels(clust_vec)) {
            enrichment_list[[length(enrichment_list)+1]] <- data.frame(Cluster=c_lvl, Feature=v, Score=(mean(vec[clust_vec==c_lvl], na.rm=T)-g_m)/g_sd)
          }
        } else {
          vec <- as.factor(vec); tbl <- table(vec, clust_vec)
          p_val <- tryCatch(chisq.test(tbl)$p.value, error=function(e) NA)
          if (!is_mandatory && (is.na(p_val) || p_val >= conf$p_cutoff)) next
          E <- rowSums(tbl) %o% colSums(tbl) / sum(tbl)
          resids <- (tbl - E) / sqrt(E)
          for (lvl in rownames(resids)) {
            for (c_lvl in colnames(resids)) {
              enrichment_list[[length(enrichment_list)+1]] <- data.frame(Cluster=c_lvl, Feature=sprintf("%s (%s)", v, lvl), Score=resids[lvl, c_lvl])
            }
          }
        }
      }

      # 3. RENDER COMPLEXHEATMAP
      if (length(enrichment_list) > 0) {
        heat_df <- bind_rows(enrichment_list)

        # --- NEW: Transpose Matrix (Clusters = Cols, Features = Rows) ---
        mat_wide <- heat_df %>%
          pivot_wider(names_from=Cluster, values_from=Score, values_fill=0) %>%
          column_to_rownames("Feature")

        safe_cat_name <- gsub("_", " ", cat_name)
        heatmap_file <- file.path(pheno_dir, paste0(job_id, "_", algo_name, "_Heatmap_", cat_name, ".png"))

        png(heatmap_file, width=1200, height=600)

        ht <- ComplexHeatmap::Heatmap(
          as.matrix(mat_wide),
          name = "Enrichment",
          column_title = sprintf("Clinical Enrichment: %s\nJob: %s | Algorithm: %s | (p < %s)",
                                 safe_cat_name, job_id, toupper(algo_name), conf$p_cutoff),

          # --- NEW: Axis and Label formatting ---
          column_names_rot = 0,         # Clusters straight on X-axis
          column_names_side = "bottom", # Clusters at the bottom
          row_names_side = "right",      # Clinical parameters on the left Y-axis
          row_names_gp = gpar(fontsize = 14),

          clustering_distance_rows = "euclidean",
          cluster_columns = FALSE,
          col = circlize::colorRamp2(c(-2, 0, 2), c("#2980b9", "white", "#c0392b"))
        )

        # --- NEW: Move the legend to the bottom ---
        ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
        dev.off()

        results$plots[[algo_name]][["heatmap"]][[cat_name]] <- ht
      }


      # --- NEW: OUTER LOOP FOR HIGHLIGHT GRIDS ---
      cluster_levels <- levels(valid_df$Cluster)
      my_comparisons <- combn(cluster_levels, 2, simplify = FALSE)

      for (hl_cat_name in names(hl_list)) {
        hl_vars <- hl_list[[hl_cat_name]]
        hl_plots <- list()

        for (hv in hl_vars) {
          if (!(hv %in% colnames(valid_df))) next

          # CONTINUOUS: Boxplots with Jitter & ggpubr brackets
          if (is.numeric(valid_df[[hv]])) {
            y_max <- max(valid_df[[hv]], na.rm = TRUE)
            y_min <- min(valid_df[[hv]], na.rm = TRUE)
            y_range <- y_max - y_min

            hl_plots[[hv]] <- ggplot(valid_df, aes(x = Cluster, y = .data[[hv]], fill = Cluster)) +
              geom_boxplot(outlier.shape = NA, alpha = 0.7) +
              geom_jitter(width = 0.2, alpha = 0.5, color = "grey30") +
              ggpubr::stat_compare_means(
                comparisons = my_comparisons,
                method = "wilcox.test",
                label = "p.format"
              ) +
              coord_cartesian(ylim = c(y_min, y_max + 0.15 * y_range)) +
              theme_project_base() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(title = hv, x = NULL)

            # CATEGORICAL: Stacked Barplots with Custom Fisher Subtitle
          } else {
            p_vals <- sapply(my_comparisons, function(p) {
              sub_df <- valid_df %>% filter(Cluster %in% p) %>% droplevels()
              tbl <- table(sub_df$Cluster, sub_df[[hv]])
              tryCatch(fisher.test(tbl, simulate.p.value = TRUE)$p.value, error = function(e) NA)
            })

            valid_idx <- !is.na(p_vals)
            fdr_vals <- rep(NA, length(p_vals))
            fdr_vals[valid_idx] <- p.adjust(p_vals[valid_idx], method = "fdr")

            sig_pairs <- c()
            for (i in seq_along(fdr_vals)) {
              if (!is.na(fdr_vals[i]) && fdr_vals[i] < 0.05) {
                n1 <- gsub("Cluster ", "C", my_comparisons[[i]][1])
                n2 <- gsub("Cluster ", "C", my_comparisons[[i]][2])
                sig_pairs <- c(sig_pairs, sprintf("%s vs %s (FDR=%.2f)", n1, n2, fdr_vals[i]))
              }
            }

            subtitle_str <- if(length(sig_pairs) == 0) "Sig Pairs (FDR < 0.05): ns" else paste("Sig:", paste(sig_pairs, collapse = " | "))

            hl_plots[[hv]] <- ggplot(valid_df, aes(x=Cluster, fill=as.factor(.data[[hv]]))) +
              geom_bar(position="stack", color = "black") +
              theme_project_base() +
              theme(axis.text.x = element_text(angle=45, hjust=1),
                    plot.subtitle = element_text(size = 9, color = "grey40")) +
              labs(title=hv, subtitle=subtitle_str, x=NULL, y="Count", fill=hv)
          }
        }

        # Stitch this specific category grid
        if (length(hl_plots) > 0) {
          safe_hl_name <- gsub("_", " ", hl_cat_name)
          master_grid <- wrap_plots(hl_plots, ncol=5) +
            plot_annotation(title=sprintf("Parameters: %s | %s (%s)", safe_hl_name, job_id, toupper(algo_name)))

          results$plots[[algo_name]][["highlights"]][[hl_cat_name]] <- master_grid
          ggsave(file.path(pheno_dir, paste0(job_id, "_", algo_name, "_HL_Grid_", hl_cat_name, ".png")), master_grid, width=9, height=4, bg="white")
        }
      }
    }
    return(results)
  }
}

# 2.5 Predictive Modeling Pipeline Wrapper (LASSO)
# Description:
#   Sets up a penalized logistic regression (LASSO) using glmnet. This wrapper
#   filters the cohorts to a binary classification task, runs cross-validation
#   to find the optimal lambda, and extracts the features with non-zero weights.
# Args:
#   omics_mat: Cleaned omics matrix.
#   clin_df: Clinical metadata dataframe.
#   lasso_conf: Target column and the two groups to classify.
#   model_id: String identifier for saving outputs.
#   base_output_dir: Path for saving results.
# Returns:
#   List containing the CV curve plot, the coefficients bar chart, and weight data.
wrap_predictive_pipeline <- function(omics_mat, clin_df, lasso_conf, model_id, base_output_dir) {

  results <- list(status = "success", plots = list(), data = list())

  # Setup Directories
  pred_base_dir <- file.path(base_output_dir, "LASSO_Models")
  pred_sub_dir <- file.path(pred_base_dir, model_id)
  if (!dir.exists(pred_sub_dir)) dir.create(pred_sub_dir, recursive = TRUE)

  # Run LASSO Engine
  lasso_res <- run_lasso_engine(omics_mat, clin_df, lasso_conf$group_col,
                                lasso_conf$target_groups, title = lasso_conf$title)

  if (is.null(lasso_res)) {
    results$status <- "failed"; results$error_msg <- "Model failed to converge or insufficient groups."; return(results)
  }

  # Package Results
  results$data$lasso <- lasso_res$data
  results$plots$lasso_cv <- lasso_res$cv_plot
  results$plots$lasso_weights <- lasso_res$coeff_plot

  # Automated Saving
  write_csv(lasso_res$data, file.path(pred_sub_dir, paste0(model_id, "_LASSO_Coefficients.csv")))

  png(file.path(pred_sub_dir, paste0(model_id, "_LASSO_CV.png")), width=7, height=6, units="in", res=300)
  plot(lasso_res$cv_plot); title(main = sprintf("%s: LASSO Tuning", lasso_conf$title), line = 2.5); dev.off()

  ggsave(file.path(pred_sub_dir, paste0(model_id, "_LASSO_Weights.png")), lasso_res$coeff_plot, width=7, height=6, dpi=300, bg="white")

  return(results)
}

# WORKERS 2: Clustering & Predictive Modeling -------------------------------

# w2.1 Shared Clustering Worker Engine
# Description:
#   The mathematical core that calculates cluster assignments on true high-dimensional
#   data. Supports Classic (K-Means, PAM, H-Clust), Density (HDBSCAN), and Graph (Louvain).
run_clustering_worker <- function(patient_matrix, k, algo_name, algo_conf, is_dist_matrix = FALSE) {
  set.seed(42)
  cluster_labels <- NULL

  # NEW LOGIC: Handle pre-computed distance matrices
  if (isTRUE(is_dist_matrix)) {
    dist_obj <- as.dist(patient_matrix)
  } else {
    dist_metric <- if(!is.null(algo_conf$distance)) algo_conf$distance else "euclidean"
    dist_obj <- dist(patient_matrix, method = dist_metric)
  }

  if (algo_name == "kmeans") {
    # Kmeans CANNOT take a dist object, it needs the raw continuous coordinates.
    if (isTRUE(is_dist_matrix)) return(NULL)
    n_start <- if(!is.null(algo_conf$nstart)) algo_conf$nstart else 25
    km_res <- kmeans(patient_matrix, centers = k, nstart = n_start)
    cluster_labels <- as.factor(paste("Cluster", km_res$cluster))

  } else if (algo_name == "pam") {
    library(cluster)
    # PAM gracefully accepts dist objects if diss = TRUE
    if (isTRUE(is_dist_matrix)) {
      pam_res <- pam(dist_obj, k = k, diss = TRUE)
    } else {
      dist_metric <- if(!is.null(algo_conf$metric)) algo_conf$metric else "euclidean"
      pam_res <- pam(patient_matrix, k = k, metric = dist_metric)
    }
    cluster_labels <- as.factor(paste("Cluster", pam_res$clustering))

  } else if (algo_name == "hclust") {
    link_method <- if(!is.null(algo_conf$method)) algo_conf$method else "ward.D2"
    hc_res <- hclust(dist_obj, method = link_method)
    cluster_labels <- as.factor(paste("Cluster", cutree(hc_res, k = k)))

    # NEW: Silently attach the raw tree object as an attribute for the wrapper to find
    attr(cluster_labels, "tree") <- hc_res

  } else if (algo_name == "louvain") {
    # Louvain relies on high-dim Euclidean space for kNN calculation via dbscan.
    # Skip if passed a raw pre-computed Gower distance.
    if (isTRUE(is_dist_matrix)) return(NULL)
    # ... (Rest of your existing Louvain logic)
  }

  if (!is.null(cluster_labels)) {
    # If it was a dist matrix, rownames are attributes of the matrix
    if (isTRUE(is_dist_matrix)) names(cluster_labels) <- attr(dist_obj, "Labels") else names(cluster_labels) <- rownames(patient_matrix)
  }

  return(cluster_labels)
}

# w2.2 Algorithm-Aware Diagnostic & Clinical Correlation Engine
run_clustering_diagnostics <- function(clustering_payload, omics_mat, clin_df, dict_df, diag_conf, cluster_conf, job_id, base_output_dir) {
  library(cluster)
  library(dbscan)
  library(igraph)

  results <- list(plots = list(), paths = list(), best_algo = NULL)

  diag_dir <- file.path(base_output_dir, "Diagnostics", job_id)
  if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

  dist_mat <- clustering_payload$data$dist_matrix
  assignments <- clustering_payload$data$assignments
  patient_mat <- t(omics_mat) # Rows = Patients, Cols = Features

  if(is.null(dist_mat) || length(assignments) == 0) return(NULL)

  global_metrics_list <- list()
  indiv_conf_list <- list()

  for (algo in names(assignments)) {
    labels <- assignments[[algo]]
    algo_conf <- cluster_conf$algorithms[[algo]]

    global_grade <- ""
    indiv_scores <- numeric(length(labels))

    # A. TOPOLOGY-SPECIFIC SCORING
    # 1. CENTROID (K-Means, PAM, H-Clust) -> Silhouette Score
    if (algo %in% c("kmeans", "pam", "hclust")) {
      numeric_clusters <- as.numeric(as.factor(labels))

      if (length(unique(numeric_clusters)) > 1) {
        sil <- silhouette(numeric_clusters, dist_mat)
        global_grade <- sprintf("Silhouette: %.3f", mean(sil[, "sil_width"]))
        indiv_scores <- sil[, "sil_width"]
      } else {
        global_grade <- "Silhouette: 0.000 (1 Cluster)"
        indiv_scores <- rep(0, length(labels))
      }
      topology_type <- "Centroid"

      # 2. DENSITY (HDBSCAN) -> Noise Ratio & Membership Probability
    } else if (algo == "hdbscan") {
      min_pts <- if(!is.null(algo_conf$minPts)) algo_conf$minPts else 5
      hd_res <- dbscan::hdbscan(patient_mat, minPts = min_pts)

      noise_pct <- (sum(hd_res$cluster == 0) / length(hd_res$cluster)) * 100
      global_grade <- sprintf("Noise: %.1f%%", noise_pct)
      indiv_scores <- hd_res$membership_prob
      topology_type <- "Density"

      # 3. GRAPH (Louvain) -> Modularity & Internal Edge Ratio
    } else if (algo == "louvain") {
      knn_k <- if(!is.null(algo_conf$k_neighbors)) algo_conf$k_neighbors else 10
      res_val <- if(!is.null(algo_conf$resolution)) algo_conf$resolution else 1.0

      knn_obj <- dbscan::kNN(patient_mat, k = knn_k)
      edges <- data.frame(
        from = rep(1:nrow(patient_mat), each = knn_k),
        to = as.vector(t(knn_obj$id)),
        weight = 1 / (1 + as.vector(t(knn_obj$dist)))
      )
      g <- igraph::simplify(igraph::graph_from_data_frame(edges, directed = FALSE), remove.multiple = TRUE, remove.loops = TRUE)
      lc <- igraph::cluster_louvain(g, resolution = res_val)

      global_grade <- sprintf("Modularity: %.3f", igraph::modularity(lc))

      # Calculate internal edge ratio for each patient
      A <- igraph::as_adjacency_matrix(g, sparse=FALSE, attr="weight")
      comm <- lc$membership
      indiv_scores <- sapply(1:igraph::vcount(g), function(i) {
        my_comm <- comm[i]
        internal_edges <- sum(A[i, comm == my_comm])
        total_edges <- sum(A[i, ])
        if(total_edges == 0) return(0) else return(internal_edges / total_edges)
      })
      topology_type <- "Graph / Network"
    }

    # B. BIOLOGICAL DISTINCTNESS (Kruskal-Wallis)
    # Exclude Noise points so they don't skew the true endotype differences
    valid_idx <- labels != "Noise"
    sig_proteins <- 0

    if (length(unique(labels[valid_idx])) >= 2) {
      p_vals <- apply(patient_mat[valid_idx, ], 2, function(x) {
        kruskal.test(x ~ labels[valid_idx])$p.value
      })
      sig_proteins <- sum(p_vals < 0.05, na.rm = TRUE)
    }

    # Compile the Global Table Row
    global_metrics_list[[length(global_metrics_list) + 1]] <- data.frame(
      Algorithm = toupper(algo),
      Topology = topology_type,
      Mathematical_Grade = global_grade,
      Smallest_Cluster_N = min(table(labels[valid_idx])),
      Biological_Distinctness = sprintf("%d Sig. Proteins", sig_proteins)
    )

    # Compile the Individual Patient Confidence Data
    indiv_conf_list[[length(indiv_conf_list) + 1]] <- data.frame(
      Subject_ID = rownames(patient_mat),
      Algorithm = toupper(algo),
      Confidence_Score = indiv_scores
    )
  }


  # C. RENDER TABLE & PLOTS
  metrics_df <- bind_rows(global_metrics_list)

  # Format table to HTML natively
  results$plots[["leaderboard"]] <- metrics_df %>%
    kableExtra::kbl(format = "html", escape = FALSE, caption = sprintf("Algorithm Performance Metrics: %s", job_id)) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE) %>%
    kableExtra::column_spec(1, bold = TRUE)

  # Individual Confidence Correlation Plots
  if (!is.null(diag_conf$correlate_vars) && length(indiv_conf_list) > 0) {
    indiv_df <- bind_rows(indiv_conf_list)
    plot_df <- clin_df %>%
      select(Subject_ID, any_of(diag_conf$correlate_vars)) %>%
      inner_join(indiv_df, by = "Subject_ID")

    for (v in diag_conf$correlate_vars) {
      if (!(v %in% colnames(plot_df))) next
      var_class <- if (v %in% dict_df$Variable) dict_df$Class[dict_df$Variable == v] else "Continuous"

      p_var <- NULL
      if (grepl("Continuous|Ordinal", var_class, ignore.case = TRUE)) {
        p_var <- ggplot(plot_df, aes(x = .data[[v]], y = Confidence_Score)) +
          geom_point(aes(color = Algorithm), alpha = 0.6) +
          geom_smooth(method = "lm", color = "black", se = FALSE) +
          facet_wrap(~Algorithm, scales = "free_y") +
          theme_project_base() + theme(legend.position = "none") +
          labs(title = sprintf("Confidence vs %s", v),
               subtitle = sprintf("Job: %s | Checking for continuous phenotypic confounders", job_id),
               x = v, y = "Patient Confidence Score (Sil / Prob / Edge Ratio)")

      } else {
        p_var <- ggplot(plot_df %>% filter(!is.na(.data[[v]])), aes(x = as.factor(.data[[v]]), y = Confidence_Score, fill = as.factor(.data[[v]]))) +
          geom_boxplot(alpha = 0.5, outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
          facet_wrap(~Algorithm, scales = "free_y") +
          theme_project_base() + theme(legend.position = "none") +
          labs(title = sprintf("Confidence Profile: %s", v),
               subtitle = sprintf("Job: %s | Checking prediction confidence across groups", job_id),
               x = v, y = "Patient Confidence Score (Sil / Prob / Edge Ratio)")
      }

      if (!is.null(p_var)) {
        safe_v <- gsub("[^A-Za-z0-9_.-]", "_", v)
        path_var <- file.path(diag_dir, paste0(job_id, "_Conf_vs_", safe_v, ".png"))
        ggsave(path_var, p_var, width = 10, height = 5, bg = "white")

        results$plots[[v]] <- p_var
        results$paths[[v]] <- path_var
      }
    }
  }

  return(results)
}

# w2.3 LASSO Engine Worker
run_lasso_engine <- function(mat, clin, target_col, target_groups, title) {
  clin_sub <- clin %>% filter(!!sym(target_col) %in% target_groups)
  x <- t(mat[, clin_sub$Subject_ID])
  y <- factor(clin_sub[[target_col]])

  set.seed(42)
  cv_fit <- cv.glmnet(x, y, family = "binomial", type.measure = "auc")

  tmp_coeffs <- coef(cv_fit, s = "lambda.min")
  df_coeffs <- data.frame(Feature = rownames(tmp_coeffs), Beta = as.numeric(tmp_coeffs)) %>%
    filter(Beta != 0, Feature != "(Intercept)") %>% arrange(desc(abs(Beta)))

  p_lasso <- ggplot(df_coeffs, aes(x = reorder(Feature, Beta), y = Beta, fill = Beta > 0)) +
    geom_bar(stat = "identity", color = "white") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = HM_Z_HIGH, "FALSE" = HM_Z_LOW),
                      labels = c("Predicts Reference", "Predicts Target")) +
    theme_project_base() +
    labs(title = sprintf("%s: Predictive Weights (LASSO)", title),
         subtitle = "Features selected by penalized regression",
         x = NULL, y = "Beta Coefficient (Model Weight)",
         fill = "Direction")

  return(list(cv_plot = cv_fit, coeff_plot = p_lasso, data = df_coeffs))
}
