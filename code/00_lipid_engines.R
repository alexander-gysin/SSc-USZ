
# 00_lipid_engines.R
# Purpose: Custom analytical engines for highly sparse mass spectrometry lipidomics.
# Architecture: WRAPPERS (orchestration/routing) at the top, WORKERS (math/plots) at the bottom.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(data.table)

# PART 1: WRAPPERS Topography ---------------------------------------------

#' Wrapper: Missingness Topography
#' Coordinates the calculation of global missingness and generates diagnostic plots.
wrap_missingness_topography <- function(lipid_mat, config, output_dir, lod_df) {

  # 1. Calculate missingness stats
  missing_stats <- worker_calc_missingness(lipid_mat)

  # 2. Generate Plots
  p_bars <- worker_plot_missingness_bars(missing_stats$df)
  p_heatmap <- worker_plot_comissingness_heatmap(lipid_mat, missing_stats$df)
  p_lod_range <- worker_plot_lod_scatter(lipid_mat, missing_stats$df, lod_df)

  # 3. Dynamic Export
  ggsave(file.path(output_dir, "topography_missingness_bars.png"), plot = p_bars,
         width = config$export$width, height = config$export$height, dpi = config$export$dpi, bg = "white")
  ggsave(
    file.path(output_dir, "topography_lod_scatter_range.png"), plot = p_lod_range,
    width = config$export$width, height = config$export$height, dpi = config$export$dpi, bg = "white"
  )

  # Note: ComplexHeatmap uses base R graphics, saving requires png() device directly if needed.

  # 4. Compile output
  return(list(
    diagnostics = list(total_lipids = nrow(missing_stats$df)),
    data = list(missing_stats = missing_stats$df),
    plots = list(
      missingness_bars = p_bars,
      co_missingness_heatmap = p_heatmap,
      lod_range = p_lod_range
    )
  ))
}

#' Wrapper: Distribution Signatures
#' Coordinates extraction of strictly measured values and plots histograms.
wrap_distribution_signatures <- function(lipid_mat, config, output_dir, project_colors_func = NULL, lod_df) {

  # Pass lod_df down to the worker
  p_hists <- worker_plot_histograms(lipid_mat, config, project_colors_func, lod_df)

  ggsave(file.path(output_dir, "topography_distribution_signatures.png"), plot = p_hists,
         width = 12, height = 8, dpi = config$export$dpi, bg = "white")

  return(list(plots = list(hists = p_hists)))
}

#' Wrapper: Phenotype Missingness
#' Evaluates missingness against clinical variables using dynamic dictionary routing.
wrap_phenotype_completeness <- function(lipid_mat, clin_df, clin_dict, config, output_dir) {

  # Ensure clean subfolder for phenotypic maps
  pheno_out_dir <- file.path(output_dir, "phenotype_completeness")
  if (!dir.exists(pheno_out_dir)) dir.create(pheno_out_dir, recursive = TRUE)

  p_pheno_list <- worker_plot_phenotype_completeness(
    mat = lipid_mat,
    clin_df = clin_df,
    clin_dict = clin_dict,
    config = config,
    out_dir = pheno_out_dir
  )

  return(list(plots = p_pheno_list))
}


# PART 2: WORKERS (Data Processing & Visualization)-----------------------------

#' Worker: Calculate Missingness
worker_calc_missingness <- function(mat) {
  # Apply across rows (Lipids)
  missing_pct <- apply(mat, 1, function(x) sum(is.na(x) | x == 0) / length(x))

  df <- data.frame(Lipid = names(missing_pct), Missingness = missing_pct, row.names = NULL) %>%
    dplyr::distinct(Lipid, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(Missingness))

  return(list(df = df))
}

#' Worker: Plot Missingness Bars
worker_plot_missingness_bars <- function(missing_df) {
  missing_df$Lipid <- factor(missing_df$Lipid, levels = missing_df$Lipid)
  missing_df$is_100 <- missing_df$Missingness == 1

  p <- ggplot(missing_df, aes(x = Lipid, y = Missingness, fill = is_100)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("FALSE" = "coral", "TRUE" = "firebrick")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.05), labels = scales::percent_format(accuracy = 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray65", linewidth = 0.6),
      legend.position = "none"
    ) +
    labs(x = "Lipid Species", y = "Missingness (%)", title = "Ranked Missingness per Lipid")

  return(p)
}

#' Worker: Plot Co-Missingness Heatmap (ComplexHeatmap)
worker_plot_comissingness_heatmap <- function(mat, missing_df) {
  binary_mat <- ifelse(is.na(mat) | mat == 0, 1, 0)
  co_mat <- binary_mat %*% t(binary_mat)
  co_mat_pct <- (co_mat / ncol(mat)) * 100

  missing_ordered <- missing_df$Missingness[match(rownames(mat), missing_df$Lipid)] * 100

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Missingness = ComplexHeatmap::anno_barplot(missing_ordered, gp = grid::gpar(fill = "gray50", col = "black"), height = grid::unit(2, "cm")),
    annotation_name_side = "left"
  )

  #col_fun <- circlize::colorRamp2(c(min(co_mat_pct), 100), c("firebrick", "navy"))
  col_fun <- circlize::colorRamp2(
    c(0, 25, 50, 75, 99.9999, 100),
    c(
      "#FDE725",
      "#5EC962",
      "#21918C",
      "#3B528B",
      "#440154",
      "firebrick"
    )
  )

  p <- ComplexHeatmap::Heatmap(
    co_mat_pct, name = "Co-Missingness\n(% Samples)", col = col_fun, top_annotation = ha,
    cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE,
    column_names_rot = 45, column_names_side = "bottom", column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 9),
    show_column_dend = FALSE, row_dend_width = grid::unit(4, "cm"), row_title = "Clustered Lipids", column_title = "Clustered Lipids"
  )
  return(p)
}

#' Worker: Limit of Detection (LOD) Scatter
worker_plot_lod_scatter <- function(mat, missing_df, lod_df) {
  lod_vec <- setNames(lod_df$LOD_Matrix,
                      lod_df$Mediator)
  mean_lod_distance <- sapply(rownames(mat), function(lipid) {
    x <- mat[lipid, ]
    measured <- x[!is.na(x) & x > 0]
    lod <- lod_vec[lipid]
    if(length(measured) > 0 && !is.na(lod)) {
      mean(log2(measured / lod))
    } else {
      NA_real_
    }
  })

  plot_long <- as.data.table(mat, keep.rownames = "Lipid")
  plot_long <- melt(
    plot_long,
    id.vars = "Lipid",
    variable.name = "Sample",
    value.name = "Abundance"
  )
  plot_long[, LOD := lod_vec[Lipid]]
  plot_long <- plot_long[
    !is.na(Abundance) &
      !is.na(LOD) &
      Abundance > 0
  ]
  plot_long[, Log2FC_LOD := log2(Abundance / LOD)]
  plot_long <- merge(
    plot_long,
    missing_df,
    by = "Lipid",
    all.x = TRUE
  )

  int_df <- data.frame(Lipid = names(mean_lod_distance), Mean_LOD_FC = mean_lod_distance)
  plot_df <- dplyr::left_join(missing_df, int_df, by = "Lipid")
  # Top 5 lipids by missingness
  label_df <- plot_df |>
    dplyr::filter(!is.na(Mean_LOD_FC), !is.na(Missingness)) |>
    dplyr::distinct(Lipid, .keep_all = TRUE)

  summary_df <- plot_long[
    ,
    .(
      min_fc = min(Log2FC_LOD),
      max_fc = max(Log2FC_LOD),
      mean_fc = mean(Log2FC_LOD)
    ),
    by = .(Lipid, Missingness)
  ]

  plot <- ggplot() +

    # 1. BACKGROUND: Individual raw measurements (from plot_long)
    geom_jitter(
      data = plot_long,
      aes(x = Log2FC_LOD, y = Missingness),
      colour = "grey70",
      alpha = 0.4,
      size = 1.2,
      height = 0.01 # Slight vertical jitter prevents points from perfectly overlapping
    ) +

    # 2. MIDGROUND: The Min-Max Range Segment (from summary_df)
    geom_segment(
      data = summary_df,
      aes(x = min_fc, xend = max_fc, y = Missingness, yend = Missingness),
      alpha = 0.8,
      colour = "grey40",
      linewidth = 0.8
    ) +

    # 3. FOREGROUND: The Mean Point (from summary_df)
    geom_point(
      data = summary_df,
      aes(x = mean_fc, y = Missingness),
      colour = "navy",
      size = 3
    ) +

    # Formatting
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_minimal() +
    labs(
      title = "Limit of Detection (LOD) Topography",
      x = "LOD Foldchange (log2) of Successfully Measured Samples",
      y = "Missingness (%)"
    ) +

    # 4. LABELS
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = Mean_LOD_FC, y = Missingness, label = Lipid),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      color = "black"
    ) +

    # Theme
    theme(
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"))

  return(plot)
}

#' Worker: Plot Histograms
worker_plot_histograms <- function(mat, config, project_colors_func, lod_df) {

  fill_color <- if(!is.null(project_colors_func)) project_colors_func("primary") else "steelblue"

  # 1. Format the main data
  df_long <- as.data.frame(t(mat)) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(-SampleID, names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::filter(!is.na(Intensity) & Intensity > 0)

  # 2. Map the LOD values to each lipid
  # Assumes lod_df has 'Mediator' (lipid name) and 'LOD_Matrix'
  lod_vec <- setNames(lod_df$LOD_Matrix, lod_df$Mediator)
  df_long$LOD <- lod_vec[df_long$Lipid]

  # 3. Log transformation handling (Applies to BOTH intensity and LOD)
  x_label <- "Intensity (Raw)"
  if(isTRUE(config$missingness$log_transform_viz)) {
    df_long$Intensity <- log2(df_long$Intensity)
    df_long$LOD <- log2(df_long$LOD) # Must log transform the LOD line too!
    x_label <- "Log2(Intensity)"
  }

  # 4. Extract unique LODs for the vline (to prevent drawing 100 overlapping lines per facet)
  lod_lines <- df_long %>%
    dplyr::distinct(Lipid, LOD) %>%
    dplyr::filter(!is.na(LOD))

  # 5. Build the Plot (n=1 samples are now included natively in the histogram)
  p <- ggplot(df_long) +
    geom_histogram(
      aes(x = Intensity),
      bins = config$missingness$hist_bins,
      fill = fill_color,
      color = "black",
      alpha = 0.7
    ) +
    geom_vline(
      data = lod_lines,
      aes(xintercept = LOD),
      color = "firebrick",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 7, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(x = x_label, y = "Count") +
    facet_wrap(~Lipid, scales = "free", ncol = 8)

  return(p)
}

#' Worker: Phenotype Completeness
#' Dynamically routes continuous vs categorical clinical features via the data dictionary.
worker_plot_phenotype_completeness <- function(mat, clin_df, clin_dict, config, out_dir) {
  plots <- list()

  df_long <- as.data.frame(t(mat)) %>%
    tibble::rownames_to_column("Subject_ID") %>%
    tidyr::pivot_longer(-Subject_ID, names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::mutate(is_missing = ifelse(is.na(Intensity) | Intensity == 0, 1, 0))

  df_join <- df_long %>% dplyr::left_join(clin_df, by = "Subject_ID")

  for(v_name in config$phenotype$vars_to_test) {
    if(!v_name %in% colnames(df_join)) next

    # 1. Dynamic Dictionary Lookup
    var_class <- clin_dict %>% dplyr::filter(Variable == v_name) %>% dplyr::pull(Class) %>% tolower()
    if(length(var_class) == 0) var_class <- "unknown"

    plot_col <- v_name

    # 2. Dynamic Continuous Binning
    if(var_class %in% c("numeric", "integer", "double", "continuous")) {
      binned_col_name <- paste0(v_name, "_binned")
      n_bins <- config$phenotype$continuous_bins

      df_join <- df_join %>%
        dplyr::mutate(
          !!binned_col_name := dplyr::ntile(!!rlang::sym(v_name), n_bins),
          !!binned_col_name := paste0("Quantile ", !!rlang::sym(binned_col_name))
        )
      plot_col <- binned_col_name
    }

    # 3. Calculate and Plot
    # EXPLICIT NA FILTER: Remove any rows where the clinical variable is missing
    df_valid <- df_join %>% dplyr::filter(!is.na(!!rlang::sym(plot_col)) & !!rlang::sym(plot_col) != "Quantile NA")

    # Get the total denominator (Total number of valid subjects in the cohort for this feature)
    total_subjects <- length(unique(df_valid$Subject_ID))

    summary_df <- df_valid %>%
      dplyr::group_by(Lipid, !!rlang::sym(plot_col)) %>%
      dplyr::summarise(
        # How many successfully measured samples exist in this specific group?
        measured_count = sum(is_missing == 0),
        # Divide by TOTAL subjects to get this group's contribution to OVERALL completeness
        comp_contrib = measured_count / total_subjects,
        .groups = "drop"
      )

    # Plotting using the new comp_contrib fraction
    p <- ggplot(summary_df, aes(x = Lipid, y = comp_contrib, fill = as.factor(!!rlang::sym(plot_col)))) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top",
        legend.title = element_blank()
      ) +
      labs(
        title = paste("Completeness Partitioned by", v_name),
        y = "Overall Completeness (%)",
        x = "Lipid Species"
      )

    # 4. Dynamic Export
    file_name <- file.path(out_dir, paste0("pheno_completeness_", v_name, ".png"))
    ggsave(filename = file_name, plot = p, width = config$export$width, height = config$export$height, dpi = config$export$dpi, bg = "white")

    plots[[v_name]] <- p
  }
  return(plots)
}

# PART 3: WRAPPERS (Statistical Routing) ---------------------------------------

#' Wrapper: Sparse Differential Expression Analysis
#' Orchestrates the binarized statistical testing, handles FDR correction,
#' and dynamically routes output creation.
wrap_sparse_dea <- function(mat_bin, clin_df, config, engine, job_name, out_dir, assets_dir, project_colors_func) {

  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_df <- clin_df %>%
      dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups))
  }

  common_samples <- intersect(colnames(mat_bin), clin_df$Subject_ID)

  if (length(common_samples) < 10) {
    return(list(status = "error", error_msg = "Fewer than 10 overlapping subjects found."))
  }

  mat_bin_sub <- mat_bin[, common_samples, drop = FALSE]

  clin_sub <- clin_df %>% filter(Subject_ID %in% common_samples)
  clin_sub <- clin_sub[match(common_samples, clin_sub$Subject_ID), ]

  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_vec_stat <- ifelse(clin_sub[[config$test_var]] %in% config$target_groups, "Target", "Reference")
  } else {
    clin_vec_stat <- clin_sub[[config$test_var]]
  }

  clin_vec_plot <- clin_sub[[config$test_var]]

  results_list <- lapply(rownames(mat_bin_sub), function(lipid) {
    lipid_vec <- mat_bin_sub[lipid, ]
    if (engine == "fisher") res <- worker_run_fishers_exact(lipid_vec, clin_vec_stat)
    else if (engine == "firth") res <- worker_run_firth_logistic(lipid_vec, clin_vec_stat)

    return(cbind(Lipid = lipid, res))
  })

  raw_results_df <- do.call(rbind, results_list)
  failed_df <- raw_results_df %>% dplyr::filter(is.na(P_Value))

  results_df <- raw_results_df %>%
    dplyr::filter(!is.na(P_Value)) %>%
    dplyr::mutate(FDR = p.adjust(P_Value, method = "BH")) %>%
    dplyr::arrange(P_Value)

  fc_thresh <- if(!is.null(config$fc_cutoff)) config$fc_cutoff else 0
  sig_df <- results_df %>% dplyr::filter(FDR < config$fdr_cutoff & abs(Log2OR) >= fc_thresh)

  p_volc <- worker_plot_volcano(results_df, config, config$title)

  p_heat <- worker_plot_dea_heatmap(
    mat = mat_bin_sub,
    clin_vec = clin_vec_plot,
    clin_var_name = config$test_var,
    results_df = results_df,
    project_colors_func = project_colors_func,
    job_title = config$title
  )

  write.csv(results_df, file.path(out_dir, paste0(job_name, "_full_results.csv")), row.names = FALSE)
  write.csv(failed_df, file.path(out_dir, paste0(job_name, "_exclusion_report.csv")), row.names = FALSE)

  if (!is.null(p_volc)) {
    ggsave(file.path(out_dir, paste0(job_name, "_volcano.png")), plot = p_volc, width = 8, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(assets_dir, paste0(job_name, "_volcano.png")), plot = p_volc, width = 8, height = 6, dpi = 300, bg = "white")
  }


  png(file.path(out_dir, paste0(job_name, "_heatmap.png")), width = 8, height = 8, units = "in", res = 300)
  ComplexHeatmap::draw(p_heat)
  invisible(dev.off())

  png(file.path(assets_dir, paste0(job_name, "_heatmap.png")), width = 8, height = 8, units = "in", res = 300)
  ComplexHeatmap::draw(p_heat)
  invisible(dev.off())

  return(list(
    status = "success",
    tables = list(full_results = results_df, significant_features = sig_df, failed_features = failed_df),
    plots = list(volcano = p_volc, heatmap = p_heat)
  ))
}

# PART 4: WORKERS  -----------------------------------------------------------

#' Worker: Fisher's Exact Test
worker_run_fishers_exact <- function(lipid_vec, clin_vec) {

  clin_vec <- droplevels(as.factor(clin_vec))

  if (length(unique(na.omit(clin_vec))) < 2 || length(unique(na.omit(lipid_vec))) < 2) {
    return(data.frame(P_Value = NA_real_, Estimate_OR = NA_real_, Log2OR = NA_real_, Method = "Fisher (No Variance)"))
  }

  tbl <- table(Lipid = lipid_vec, Pheno = clin_vec)
  ft <- fisher.test(tbl)

  p_val <- ft$p.value

  if (!is.null(ft$estimate)) {
    or_val <- as.numeric(ft$estimate)
    log2_or <- log2(or_val + 1e-9)
  } else {
    or_val <- NA_real_
    log2_or <- NA_real_
  }

  return(data.frame(
    P_Value = p_val,
    Estimate_OR = or_val,
    Log2OR = log2_or,
    Method = "Fisher"
  ))
}

#' Worker: Firth's Penalized Logistic Regression
worker_run_firth_logistic <- function(lipid_vec, clin_vec) {

  df <- data.frame(Y = as.numeric(lipid_vec), X = as.numeric(clin_vec)) %>% na.omit()

  if(nrow(df) < 10 || length(unique(df$Y)) < 2) {
    return(data.frame(P_Value = NA_real_, Estimate_OR = NA_real_, Log2OR = NA_real_, Method = "Firth (Low N/Var)"))
  }

  fit <- tryCatch({
    logistf::logistf(Y ~ X, data = df)
  }, error = function(e) return(NULL))

  if (is.null(fit)) {
    return(data.frame(P_Value = NA_real_, Estimate_OR = NA_real_, Log2OR = NA_real_, Method = "Firth (Convergence Fail)"))
  }

  p_val <- unname(fit$prob["X"])
  coef_x <- unname(fit$coef["X"])
  or_x <- exp(coef_x)
  log2_or <- coef_x / log(2)

  return(data.frame(
    P_Value = p_val,
    Estimate_OR = or_x,
    Log2OR = log2_or,
    Method = "Firth"
  ))
}

#' Worker: Standard Volcano Plot
worker_plot_volcano <- function(res_df, config, job_title) {

  fc_thresh  <- if(!is.null(config$fc_cutoff)) config$fc_cutoff else 1.0
  fdr_thresh <- if(!is.null(config$fdr_cutoff)) config$fdr_cutoff else 0.05

  plot_df <- res_df %>%
    dplyr::mutate(
      Plot_Log2OR = ifelse(Log2OR > 10, 10, ifelse(Log2OR < -10, -10, Log2OR)),
      NegLog10FDR = -log10(FDR + 1e-300),
      Is_Significant = ifelse(FDR < fdr_thresh & abs(Log2OR) >= fc_thresh, "Significant", "Not Significant")
    )

  top_hits <- plot_df %>%
    dplyr::filter(Is_Significant == "Significant") %>%
    dplyr::slice_max(order_by = NegLog10FDR, n = 10)

  p <- ggplot(plot_df, aes(x = Plot_Log2OR, y = NegLog10FDR)) +
    geom_point(data = subset(plot_df, Is_Significant == "Not Significant"), color = "grey70", alpha = 0.5, size = 1.5) +
    geom_point(data = subset(plot_df, Is_Significant == "Significant"), color = "firebrick", alpha = 0.8, size = 2.5) +

    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "navy", alpha = 0.6) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "navy", alpha = 0.6) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +

    ggrepel::geom_text_repel(data = top_hits, aes(label = Lipid), size = 3.5, box.padding = 0.5, max.overlaps = Inf) +

    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14), axis.title = element_text(size = 12), panel.grid.minor = element_blank()) +

    # Dynamic Title Injection
    labs(title = sprintf("Volcano (FDR): %s", job_title), x = "Log2(Odds Ratio)", y = "-Log10(FDR)")

  return(p)
}

#' Worker: Sparse DEA Heatmap
worker_plot_dea_heatmap <- function(mat, clin_vec, clin_var_name, results_df, project_colors_func, job_title) {

  ranked_lipids <- results_df %>% dplyr::filter(!is.na(P_Value)) %>% dplyr::arrange(P_Value) %>% dplyr::pull(Lipid)
  ranked_lipids <- intersect(ranked_lipids, rownames(mat))
  plot_mat <- mat[ranked_lipids, , drop = FALSE]

  if (is.numeric(clin_vec)) {
    low_col <- if(exists("COLOR_AGE_LOW")) COLOR_AGE_LOW else "white"
    high_col <- if(exists("COLOR_AGE_HIGH")) COLOR_AGE_HIGH else "purple4"
    col_map <- circlize::colorRamp2(
      c(min(clin_vec, na.rm = TRUE), max(clin_vec, na.rm = TRUE)),
      c(low_col, high_col)
    )
  } else {
    unique_levels <- as.character(na.omit(unique(clin_vec)))
    col_map <- project_colors_func(unique_levels)
  }

  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(Pheno = clin_vec),
    col = list(Pheno = col_map),
    na_col = "grey80",
    annotation_name_side = "left",
    annotation_legend_param = list(Pheno = list(title = clin_var_name))
  )

  col_fun <- circlize::colorRamp2(c(0, 1), c("grey90", "firebrick"))

  ht <- ComplexHeatmap::Heatmap(
    plot_mat,
    name = "Detection",
    col = col_fun,
    top_annotation = ha_col,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    row_title = sprintf("Tested Lipids (n=%d)", length(ranked_lipids)),

    # Dynamic Title Injection (adds Job Title to the top of the Heatmap footnote)
    column_title = sprintf("%s\nGrouped by %s\n\n*Note: Lipids with zero variance across groups are excluded.", job_title, clin_var_name),
    column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    border = TRUE
  )

  return(ht)
}



# PART 5: Raw Abundances --------------------------------------------------

#' Wrapper: Raw Abundance Plotting
#' Merges contrast groups, calculates conditional Wilcoxon FDR, and triggers 4-col patchwork.
wrap_raw_abundances <- function(mat_raw, clin_df, config, job_name, out_dir, assets_dir, project_colors_func) {

  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_df <- clin_df %>% dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups))
  }

  common_samples <- intersect(colnames(mat_raw), clin_df$Subject_ID)
  if (length(common_samples) < 5) return(list(status = "error", error_msg = "Insufficient subjects."))

  mat_sub <- mat_raw[, common_samples, drop = FALSE]
  clin_sub <- clin_df %>% filter(Subject_ID %in% common_samples)
  clin_sub <- clin_sub[match(common_samples, clin_sub$Subject_ID), ]

  is_cont <- is.numeric(clin_sub[[config$test_var]])

  # 1. Group Merging & Formatting
  if (!is_cont) {
    target_name <- paste(config$target_groups, collapse = "+")
    ref_name <- paste(config$ref_groups, collapse = "+")

    clin_sub$Plot_Group <- dplyr::case_when(
      clin_sub[[config$test_var]] %in% config$target_groups ~ target_name,
      clin_sub[[config$test_var]] %in% config$ref_groups ~ ref_name,
      TRUE ~ NA_character_
    )
    clin_vec_plot <- clin_sub$Plot_Group
  } else {
    target_name <- "High"
    ref_name <- "Low"
    clin_vec_plot <- clin_sub[[config$test_var]]
  }

  # 2. Statistical Testing (Conditional Abundance)
  stats_list <- lapply(rownames(mat_sub), function(lipid) {
    val_vec <- mat_sub[lipid, ]

    # Isolate strictly detected data
    valid_idx <- !is.na(val_vec) & val_vec > 0 & !is.na(clin_vec_plot)
    v_valid <- val_vec[valid_idx]
    g_valid <- clin_vec_plot[valid_idx]

    if (!is_cont) {
      # Mann-Whitney U for Categorical Contasts
      v_target <- v_valid[g_valid == target_name]
      v_ref <- v_valid[g_valid == ref_name]
      if (length(v_target) >= 3 && length(v_ref) >= 3) {
        p_val <- suppressWarnings(wilcox.test(v_target, v_ref, exact = FALSE)$p.value)
      } else { p_val <- NA_real_ }
    } else {
      # Spearman Correlation for Continuous Contasts
      if (length(v_valid) >= 5) {
        p_val <- suppressWarnings(cor.test(v_valid, g_valid, method = "spearman")$p.value)
      } else { p_val <- NA_real_ }
    }
    return(data.frame(Lipid = lipid, P_Value = p_val, stringsAsFactors = FALSE))
  })

  stats_df <- do.call(rbind, stats_list)
  stats_df$FDR <- p.adjust(stats_df$P_Value, method = "BH")

  # 3. Call Plotting Worker
  p_raw <- worker_plot_raw_intensities(
    mat_raw = mat_sub,
    clin_vec = clin_vec_plot,
    stats_df = stats_df,
    is_cont = is_cont,
    target_name = target_name,
    ref_name = ref_name,
    config = config,
    project_colors_func = project_colors_func
  )

  # 4. Save dynamically scaled image (4 columns wide)
  if (!is.null(p_raw)) {
    grid_height <- max(3, ceiling(nrow(mat_sub) / 3) * 3.5)
    file_name <- paste0(job_name, "_raw_intensities.png")
    ggsave(file.path(out_dir, file_name), plot = p_raw, width = 16, height = grid_height, dpi = 300, bg = "white")
    ggsave(file.path(assets_dir, file_name), plot = p_raw, width = 16, height = grid_height, dpi = 300, bg = "white")
  }

  return(list(status = "success", plot = p_raw))
}

#' Worker: Plot Raw Intensities (Patchwork Grid)
worker_plot_raw_intensities <- function(mat_raw, clin_vec, stats_df, is_cont, target_name, ref_name, config, project_colors_func) {

  target_lipids <- sort(rownames(mat_raw))
  if (length(target_lipids) == 0) return(NULL)

  plot_list <- list()

  if (is_cont) {
    plot_x <- as.factor(dplyr::ntile(clin_vec, 3))
    levels(plot_x) <- c("Low", "Mid", "High")
  } else {
    plot_x <- factor(clin_vec, levels = c(ref_name, target_name))

    # Pass both groups simultaneously to properly coordinate the fallback palette
    base_levels <- c(config$target_groups[1], config$ref_groups[1])
    assigned_colors <- project_colors_func(base_levels)
    col_map <- setNames(assigned_colors, c(target_name, ref_name))
  }

  for (lipid in target_lipids) {
    df <- data.frame(Pheno = plot_x, Value = mat_raw[lipid, ])
    df_clean <- df %>% dplyr::filter(!is.na(Value) & Value > 0)

    if (nrow(df_clean) == 0) next

    df_clean$Log10Value <- log10(df_clean$Value)

    # Dynamic Headers (FDR + N counts)
    fdr_val <- stats_df$FDR[stats_df$Lipid == lipid]
    test_str <- if(is_cont) "Spearman" else "MWU"
    fdr_str <- if(is.na(fdr_val)) sprintf("%s FDR: NA", test_str) else sprintf("%s FDR: %.2f", test_str, fdr_val)

    if (!is_cont) {
      g_counts <- table(df_clean$Pheno)
      n_t <- if (target_name %in% names(g_counts)) g_counts[[target_name]] else 0
      n_r <- if (ref_name %in% names(g_counts)) g_counts[[ref_name]] else 0
      sub_str <- sprintf("%s | %s (n=%d) vs %s (n=%d)", fdr_str, target_name, n_t, ref_name, n_r)
      p_cols <- col_map
    } else {
      sub_str <- sprintf("%s | Total (n=%d)", fdr_str, nrow(df_clean))
      p_cols <- setNames(viridis::viridis(3), c("Low", "Mid", "High"))
    }

    p <- ggplot(df_clean, aes(x = Pheno, y = Log10Value, color = Pheno)) +
      geom_jitter(width = 0.2, alpha = 0.8, size = 2) +
      scale_color_manual(values = p_cols) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 12),
        panel.grid.minor = element_blank()
      ) +
      labs(title = lipid, subtitle = sub_str, y = "Log10(Intensity)")

    plot_list[[lipid]] <- p
  }

  if (length(plot_list) == 0) return(NULL)

  # 4 Column Patchwork Layout
  grid_plot <- patchwork::wrap_plots(plot_list, ncol = 3) +
    patchwork::plot_annotation(
      title = sprintf("Raw Conditional Abundances: %s", config$title),
      subtitle = "Values log10-transformed after excluding non-detected features.",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  return(grid_plot)
}

# PART 6: Lipid Set Enrichment Analysis ----------------------------------------

#' Wrapper: Lipid Set Enrichment Analysis (Binary Burden)
#' Evaluates the sum of detected lipids per structural family using Wilcoxon/Spearman.
wrap_lsea <- function(burden_mat, clin_df, config, job_name, out_dir, assets_dir, project_colors_func) {

  # 1. Filter to requested contrast groups
  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_df <- clin_df %>% dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups))
  }

  common_samples <- intersect(colnames(burden_mat), clin_df$Subject_ID)
  if (length(common_samples) < 5) return(list(status = "error", error_msg = "Insufficient subjects for LSEA."))

  mat_sub <- burden_mat[, common_samples, drop = FALSE]
  clin_sub <- clin_df %>% filter(Subject_ID %in% common_samples)
  clin_sub <- clin_sub[match(common_samples, clin_sub$Subject_ID), ]

  is_cont <- is.numeric(clin_sub[[config$test_var]])

  # 2. Group Merging & Formatting
  if (!is_cont) {
    target_name <- paste(config$target_groups, collapse = "+")
    ref_name <- paste(config$ref_groups, collapse = "+")

    clin_sub$Plot_Group <- dplyr::case_when(
      clin_sub[[config$test_var]] %in% config$target_groups ~ target_name,
      clin_sub[[config$test_var]] %in% config$ref_groups ~ ref_name,
      TRUE ~ NA_character_
    )
    clin_vec_plot <- clin_sub$Plot_Group
  } else {
    target_name <- "High"
    ref_name <- "Low"
    clin_vec_plot <- clin_sub[[config$test_var]]
  }

  # 3. Statistical Testing (Binary Burden)
  stats_list <- lapply(rownames(mat_sub), function(family) {
    v_valid <- mat_sub[family, ]
    g_valid <- clin_vec_plot

    # Skip families with zero variance (e.g., no lipids detected in either group)
    if (length(unique(v_valid)) <= 1) {
      return(data.frame(Family = family, P_Value = NA_real_, stringsAsFactors = FALSE))
    }

    if (!is_cont) {
      v_target <- v_valid[g_valid == target_name]
      v_ref <- v_valid[g_valid == ref_name]
      if (length(v_target) >= 3 && length(v_ref) >= 3) {
        p_val <- suppressWarnings(wilcox.test(v_target, v_ref, exact = FALSE)$p.value)
      } else { p_val <- NA_real_ }
    } else {
      if (length(v_valid) >= 5) {
        p_val <- suppressWarnings(cor.test(v_valid, g_valid, method = "spearman")$p.value)
      } else { p_val <- NA_real_ }
    }
    return(data.frame(Family = family, P_Value = p_val, stringsAsFactors = FALSE))
  })

  stats_df <- do.call(rbind, stats_list)
  stats_df$FDR <- p.adjust(stats_df$P_Value, method = "BH")

  # 4. Save Statistical Results
  write.csv(stats_df, file.path(out_dir, paste0(job_name, "_LSEA_results.csv")), row.names = FALSE)

  # 5. Call Plotting Worker
  p_lsea <- worker_plot_lsea(
    burden_mat = mat_sub,
    clin_vec = clin_vec_plot,
    stats_df = stats_df,
    is_cont = is_cont,
    target_name = target_name,
    ref_name = ref_name,
    config = config,
    project_colors_func = project_colors_func
  )

  # 6. Save dynamically scaled image (4 columns wide)
  if (!is.null(p_lsea)) {
    grid_height <- max(4, ceiling(nrow(mat_sub) / 4) * 4)
    file_name <- paste0(job_name, "_LSEA_burden.png")
    ggsave(file.path(out_dir, file_name), plot = p_lsea, width = 16, height = grid_height, dpi = 300, bg = "white")
    ggsave(file.path(assets_dir, file_name), plot = p_lsea, width = 16, height = grid_height, dpi = 300, bg = "white")
  }

  return(list(status = "success", plot = p_lsea, tables = list(lsea_stats = stats_df)))
}

#' Worker: Plot LSEA (Patchwork Boxplots + Jitter)
worker_plot_lsea <- function(burden_mat, clin_vec, stats_df, is_cont, target_name, ref_name, config, project_colors_func) {

  target_families <- sort(rownames(burden_mat))
  if (length(target_families) == 0) return(NULL)

  plot_list <- list()

  if (is_cont) {
    plot_x <- as.factor(dplyr::ntile(clin_vec, 3))
    levels(plot_x) <- c("Low", "Mid", "High")
  } else {
    plot_x <- factor(clin_vec, levels = c(ref_name, target_name))

    # Safe Color Mapping: Inherit color from the first group in the vector
    base_levels <- c(config$target_groups[1], config$ref_groups[1])
    assigned_colors <- project_colors_func(base_levels)
    col_map <- setNames(assigned_colors, c(target_name, ref_name))
  }

  for (family in target_families) {
    df <- data.frame(Pheno = plot_x, Burden = burden_mat[family, ])
    df_clean <- df %>% dplyr::filter(!is.na(Burden) & !is.na(Pheno))

    if (nrow(df_clean) == 0) next

    # Dynamic Headers (FDR + N counts)
    fdr_val <- stats_df$FDR[stats_df$Family == family]
    test_str <- if(is_cont) "Spearman" else "MWU"
    fdr_str <- if(is.na(fdr_val)) sprintf("%s FDR: NA", test_str) else sprintf("%s FDR: %.2f", test_str, fdr_val)

    if (!is_cont) {
      g_counts <- table(df_clean$Pheno)
      n_t <- if (target_name %in% names(g_counts)) g_counts[[target_name]] else 0
      n_r <- if (ref_name %in% names(g_counts)) g_counts[[ref_name]] else 0
      sub_str <- sprintf("%s | %s (n=%d) vs %s (n=%d)", fdr_str, target_name, n_t, ref_name, n_r)
      p_cols <- col_map
    } else {
      sub_str <- sprintf("%s | Total (n=%d)", fdr_str, nrow(df_clean))
      p_cols <- setNames(viridis::viridis(3), c("Low", "Mid", "High"))
    }

    # Boxplot + Jitter visualization
    p <- ggplot(df_clean, aes(x = Pheno, y = Burden, fill = Pheno, color = Pheno)) +
      geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
      geom_jitter(width = 0.18, height = 0.05, alpha = 0.8, size = 2) +
      scale_fill_manual(values = p_cols) +
      scale_color_manual(values = p_cols) +
      scale_y_continuous(breaks = function(limits) seq(0, ceiling(limits[2]), by = 1)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "grey30"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 10),
        panel.grid.minor = element_blank()
      ) +
      labs(title = family, subtitle = sub_str, y = "Detected Lipids (Count)")

    plot_list[[family]] <- p
  }

  if (length(plot_list) == 0) return(NULL)

  # 4 Column Patchwork Layout
  grid_plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = sprintf("LSEA Binary Burden: %s", config$title),
      subtitle = "Boxplots representing the total count of detected lipids per biological family.",
      theme = theme(plot.title = element_text(size = 17, face = "bold"))
    )

  return(grid_plot)
}
