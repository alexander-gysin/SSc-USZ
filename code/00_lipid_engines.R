
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
wrap_sparse_dea <- function(mat, clin_df, config, engine, job_name, out_dir) {

  # --- NEW: Explicitly isolate the requested contrast groups ---
  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_df <- clin_df %>%
      dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups))
  }

  # 1. Subject Alignment & Validation
  common_samples <- intersect(colnames(mat), clin_df$Subject_ID)

  if (length(common_samples) < 10) {
    return(list(status = "error", error_msg = "Fewer than 10 overlapping subjects found."))
  }

  mat_sub <- mat[, common_samples, drop = FALSE]
  clin_sub <- clin_df %>% filter(Subject_ID %in% common_samples)

  # Ensure exact ordering to prevent metadata mismatch
  clin_sub <- clin_sub[match(common_samples, clin_sub$Subject_ID), ]
  clin_vec <- clin_sub[[config$test_var]]

  # 2. Iterate and Execute Statistical Workers
  results_list <- lapply(rownames(mat_sub), function(lipid) {
    lipid_vec <- mat_sub[lipid, ]

    # Route to the correct mathematical engine
    if (engine == "fisher") {
      res <- worker_run_fishers_exact(lipid_vec, clin_vec)
    } else if (engine == "firth") {
      res <- worker_run_firth_logistic(lipid_vec, clin_vec)
    } else {
      stop("Unknown statistical engine requested.")
    }

    # Prepend Lipid name
    res <- cbind(Lipid = lipid, res)
    return(res)
  })

  # 3. Post-Processing & Multiple Testing Correction
  results_df <- do.call(rbind, results_list) %>%
    # Drop lipids that failed mathematically (e.g., zero variance in this specific slice)
    dplyr::filter(!is.na(P_Value)) %>%
    dplyr::mutate(
      FDR = p.adjust(P_Value, method = "BH")
    ) %>%
    dplyr::arrange(P_Value)

  sig_df <- results_df %>% dplyr::filter(FDR < config$p_cutoff)

  # 4. Visualization Routing
  p_volc <- worker_plot_volcano(results_df, config, job_name)

  # 5. Dynamic Export
  write.csv(results_df, file.path(out_dir, paste0(job_name, "_full_results.csv")), row.names = FALSE)

  if (!is.null(p_volc)) {
    ggsave(
      filename = file.path(out_dir, paste0(job_name, "_volcano.png")),
      plot = p_volc,
      width = 8, height = 6, dpi = 300, bg = "white"
    )
  }

  # 6. Compile Output
  return(list(
    status = "success",
    tables = list(
      full_results = results_df,
      significant_features = sig_df
    ),
    plots = list(volcano = p_volc)
  ))
}


# PART 4: WORKERS (Statistical Engines & DEA Plots) ----------------------------

#' Worker: Fisher's Exact Test
#' Calculates probability for categorical clinical variables vs binarized lipids.
worker_run_fishers_exact <- function(lipid_vec, clin_vec) {

  # Ensure factors drop unrepresented levels to strictly enforce 2x2 if filtered
  clin_vec <- droplevels(as.factor(clin_vec))

  if (length(unique(na.omit(clin_vec))) < 2 || length(unique(na.omit(lipid_vec))) < 2) {
    return(data.frame(P_Value = NA_real_, Estimate_OR = NA_real_, Log2OR = NA_real_, Method = "Fisher (No Variance)"))
  }

  tbl <- table(Lipid = lipid_vec, Pheno = clin_vec)
  ft <- fisher.test(tbl)

  p_val <- ft$p.value

  # --- NEW: Safely extract odds ratio ONLY if it exists (2x2 tables) ---
  if (!is.null(ft$estimate)) {
    or_val <- as.numeric(ft$estimate)
    log2_or <- log2(or_val + 1e-9)
  } else {
    or_val <- NA_real_     # Fallback for >2x2 tables
    log2_or <- NA_real_    # Fallback for >2x2 tables
  }

  return(data.frame(
    P_Value = p_val,
    Estimate_OR = or_val,
    Log2OR = log2_or,
    Method = "Fisher"
  ))
}

#' Worker: Firth's Penalized Logistic Regression
#' Calculates probability for continuous clinical variables, penalizing likelihood to prevent
#' infinite standard errors caused by quasi-separation in sparse features.
worker_run_firth_logistic <- function(lipid_vec, clin_vec) {

  # Combine and clean NAs for modeling
  df <- data.frame(Y = as.numeric(lipid_vec), X = as.numeric(clin_vec)) %>% na.omit()

  # Safeguard: Ensure sufficient N and variance for a likelihood model
  if(nrow(df) < 10 || length(unique(df$Y)) < 2) {
    return(data.frame(P_Value = NA_real_, Estimate_OR = NA_real_, Log2OR = NA_real_, Method = "Firth (Low N/Var)"))
  }

  # Fit Firth's Model inside tryCatch to prevent pipeline crashes on failure to converge
  fit <- tryCatch({
    logistf::logistf(Y ~ X, data = df)
  }, error = function(e) return(NULL))

  if (is.null(fit)) {
    return(data.frame(P_Value = NA_real_, Estimate_OR = NA_real_, Log2OR = NA_real_, Method = "Firth (Convergence Fail)"))
  }

  # Extract P-value (from Profile Likelihood by default in logistf)
  p_val <- unname(fit$prob["X"])

  # Extract Coefficient (which is Natural Log Odds)
  coef_x <- unname(fit$coef["X"])
  or_x <- exp(coef_x)

  # Convert Natural Log Odds to base-2 Log Odds Ratio for standard volcano comparability
  log2_or <- coef_x / log(2)

  return(data.frame(
    P_Value = p_val,
    Estimate_OR = or_x,
    Log2OR = log2_or,
    Method = "Firth"
  ))
}

#' Worker: Standard Volcano Plot
worker_plot_volcano <- function(res_df, config, job_name) {

  # 1. Clean data for plotting
  plot_df <- res_df %>%
    dplyr::mutate(
      # Cap extreme Log2OR values (e.g. infinite separation in Fisher's) for visual scale
      Plot_Log2OR = ifelse(Log2OR > 10, 10, ifelse(Log2OR < -10, -10, Log2OR)),
      NegLog10P = -log10(P_Value),
      Is_Significant = ifelse(FDR < config$p_cutoff, "Significant", "Not Significant")
    )

  # 2. Extract Top Hits for Labeling
  top_hits <- plot_df %>%
    dplyr::filter(Is_Significant == "Significant") %>%
    dplyr::slice_min(order_by = P_Value, n = 10) # Top 10 most significant

  # 3. Build Plot
  p <- ggplot(plot_df, aes(x = Plot_Log2OR, y = NegLog10P)) +

    # Background Points (Non-Significant)
    geom_point(
      data = subset(plot_df, Is_Significant == "Not Significant"),
      color = "grey70", alpha = 0.5, size = 1.5
    ) +

    # Foreground Points (Significant)
    geom_point(
      data = subset(plot_df, Is_Significant == "Significant"),
      color = "firebrick", alpha = 0.8, size = 2.5
    ) +

    # Threshold Lines
    geom_hline(
      yintercept = -log10(max(plot_df$P_Value[plot_df$FDR < config$p_cutoff], na.rm=T)),
      linetype = "dashed", color = "navy", alpha = 0.6
    ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +

    # Labels
    ggrepel::geom_text_repel(
      data = top_hits,
      aes(label = Lipid),
      size = 3.5,
      box.padding = 0.5,
      max.overlaps = Inf
    ) +

    # Formatting
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Differential Expression:", gsub("_", " ", job_name)),
      x = "Log2(Odds Ratio)",
      y = "-Log10(P-Value)"
    )

  return(p)
}
