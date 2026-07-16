# 00_lipid_engines.R
# Purpose: Custom analytical engines for highly sparse mass spectrometry lipidomics.
# Architecture: WRAPPERS (orchestration/routing) at the top, WORKERS (math/plots) at the bottom.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(ggrepel)
library(scales)

# PART 1: WRAPPERS -------------------------------------------------------------

# Wrapper: Missingness Topography
#' Coordinates the calculation of global missingness and generates diagnostic plots.
#' @param lipid_mat A numeric matrix of unimputed lipidomics (Lipids as rows, Patients as columns).
#' @param config List of visualization parameters.
#' @param output_dir Path to the output directory to export PNGs.
#' @param lod_df Data frame containing the aggregated mean LOD values per lipid.
wrap_missingness_topography <- function(lipid_mat, config, output_dir, lod_df) {

  # 1. Calculate missingness stats
  missing_stats <- worker_calc_missingness(lipid_mat)

  # 2. Generate Plots
  p_bars <- worker_plot_missingness_bars(missing_stats$df)
  p_heatmap <- worker_plot_comissingness_heatmap(lipid_mat, missing_stats$df)
  p_lod <- worker_plot_lod_scatter(lipid_mat, missing_stats$df, lod_df)

  # 3. Dynamic Export
  ggplot2::ggsave(
    filename = file.path(output_dir, "topography_missingness_bars.png"),
    plot = p_bars,
    width = config$export$width,
    height = config$export$height,
    dpi = config$export$dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(output_dir, "topography_lod_scatter.png"),
    plot = p_lod,
    width = config$export$width,
    height = config$export$height,
    dpi = config$export$dpi,
    bg = "white"
  )

  # 4. Compile output
  return(list(
    diagnostics = list(total_lipids = nrow(missing_stats$df)),
    data = list(missing_stats = missing_stats$df),
    plots = list(
      missingness_bars = p_bars,
      co_missingness_heatmap = p_heatmap,
      lod_scatter = p_lod
    )
  ))
}


# Wrapper: Distribution Signatures
#' Coordinates extraction of strictly measured values and plots histograms with LOD lines.
#' @param lipid_mat A numeric matrix of unimputed lipidomics.
#' @param config List of visualization parameters.
#' @param output_dir Path to the output directory to export PNGs.
#' @param project_colors_func Optional function to fetch project-specific colors.
#' @param lod_df Data frame containing the aggregated mean LOD values per lipid.
wrap_distribution_signatures <- function(lipid_mat, config, output_dir, project_colors_func = NULL, lod_df) {

  # Generate Plot
  p_hists <- worker_plot_histograms(lipid_mat, config, project_colors_func, lod_df)

  # Dynamic Export
  ggplot2::ggsave(
    filename = file.path(output_dir, "topography_distribution_signatures.png"),
    plot = p_hists,
    width = 12,
    height = 8,
    dpi = config$export$dpi,
    bg = "white"
  )

  return(list(
    plots = list(hists = p_hists)
  ))
}


# Wrapper: Phenotype Completeness
#' Evaluates completeness against clinical variables using dynamic dictionary routing.
#' @param lipid_mat A numeric matrix of unimputed lipidomics.
#' @param clin_df Clinical spine data frame containing subjects.
#' @param clin_dict Clinical data dictionary for variable class lookup.
#' @param config List of visualization parameters.
#' @param output_dir Path to the output directory to export PNGs.
wrap_phenotype_completeness <- function(lipid_mat, clin_df, clin_dict, config, output_dir) {

  # Ensure clean subfolder for phenotypic maps
  pheno_out_dir <- file.path(output_dir, "phenotype_completeness")
  if (!dir.exists(pheno_out_dir)) dir.create(pheno_out_dir, recursive = TRUE)

  # Generate Plots
  p_pheno_list <- worker_plot_phenotype_completeness(
    mat = lipid_mat,
    clin_df = clin_df,
    clin_dict = clin_dict,
    config = config,
    out_dir = pheno_out_dir
  )

  return(list(
    plots = p_pheno_list
  ))
}


# PART 2: WORKERS (Data Processing & Visualization)----------------------------

# Worker: Calculate Missingness
#' Safely handles both NA and 0 as missing values on row-wise matrix inputs.
worker_calc_missingness <- function(mat) {
  # Apply across rows (Lipids)
  missing_pct <- apply(mat, 1, function(x) sum(is.na(x) | x == 0) / length(x))

  df <- data.frame(
    Lipid = names(missing_pct),
    Missingness = missing_pct,
    row.names = NULL
  ) %>%
    dplyr::distinct(Lipid, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(Missingness))

  return(list(df = df))
}


# Worker: Plot Missingness Bars
#' Visualizes ranked missingness, coloring 100% missing values uniquely.
worker_plot_missingness_bars <- function(missing_df) {
  missing_df$Lipid <- factor(missing_df$Lipid, levels = missing_df$Lipid)
  missing_df$is_100 <- missing_df$Missingness == 1

  p <- ggplot2::ggplot(missing_df, ggplot2::aes(x = Lipid, y = Missingness, fill = is_100)) +
    ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "coral", "TRUE" = "firebrick")) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.05),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray80", linewidth = 0.5),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = "Lipid Species",
      y = "Missingness (%)",
      title = "Ranked Missingness per Lipid"
    )

  return(p)
}


# Worker: Plot Co-Missingness Heatmap
#' Uses ComplexHeatmap to map paired missingness combinations with hierarchical clustering.
worker_plot_comissingness_heatmap <- function(mat, missing_df) {
  binary_mat <- ifelse(is.na(mat) | mat == 0, 1, 0)
  co_mat <- binary_mat %*% t(binary_mat)
  co_mat_pct <- (co_mat / ncol(mat)) * 100

  missing_ordered <- missing_df$Missingness[match(rownames(mat), missing_df$Lipid)] * 100

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Missingness = ComplexHeatmap::anno_barplot(
      missing_ordered,
      gp = grid::gpar(fill = "gray50", col = "black"),
      height = grid::unit(2, "cm")
    ),
    annotation_name_side = "left"
  )

  # Inverted color scale: Blue for 100% (high co-missingness), Red for 0% (low)
  col_fun <- circlize::colorRamp2(c(min(co_mat_pct), 100), c("firebrick", "navy"))

  p <- ComplexHeatmap::Heatmap(
    co_mat_pct,
    name = "Co-Missingness\n(% Samples)",
    col = col_fun,
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 6),
    show_column_dend = FALSE,
    row_dend_width = grid::unit(4, "cm"),
    row_title = "Clustered Lipids",
    column_title = "Clustered Lipids"
  )

  return(p)
}


# Worker: Limit of Detection (LOD) Scatter
#' Maps LOD fold changes including underlying raw point jittering and mean metrics.
worker_plot_lod_scatter <- function(mat, missing_df, lod_df) {
  lod_vec <- setNames(lod_df$LOD_Matrix, lod_df$Mediator)

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

  # Identify top 5 lipids by Mean LOD FC to label
  label_df <- plot_df |>
    dplyr::filter(!is.na(Mean_LOD_FC)) |>
    dplyr::slice_max(Mean_LOD_FC, n = 5)

  summary_df <- plot_long[
    ,
    .(
      min_fc = min(Log2FC_LOD),
      max_fc = max(Log2FC_LOD),
      mean_fc = mean(Log2FC_LOD)
    ),
    by = .(Lipid, Missingness)
  ]

  # Multi-layer plot setup using strict local aesthetics to avoid inheritance crashes
  plot <- ggplot2::ggplot() +
    # Layer 1: Raw individual measurements
    ggplot2::geom_jitter(
      data = plot_long,
      ggplot2::aes(x = Log2FC_LOD, y = Missingness),
      colour = "grey70",
      alpha = 0.3,
      size = 1.2,
      height = 0.01
    ) +
    # Layer 2: Min-Max segments
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(x = min_fc, xend = max_fc, y = Missingness, yend = Missingness),
      alpha = 0.4,
      colour = "grey50"
    ) +
    # Layer 3: Means
    ggplot2::geom_point(
      data = summary_df,
      ggplot2::aes(x = mean_fc, y = Missingness),
      colour = "navy",
      size = 2
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format()
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Limit of Detection (LOD) Topography",
      x = "LOD Foldchange (log2) of Successfully Measured Samples",
      y = "Missingness (%)"
    ) +
    # Layer 4: Explicitly mapped text labels
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(x = Mean_LOD_FC, y = Missingness, label = Lipid),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3
    )

  return(plot)
}


# Worker: Plot Histograms
#' Plots raw or log2 facet-wrapped histograms including dashed red LOD reference thresholds.
worker_plot_histograms <- function(mat, config, project_colors_func, lod_df) {
  fill_color <- if(!is.null(project_colors_func)) project_colors_func("primary") else "steelblue"

  df_long <- as.data.frame(t(mat)) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(-SampleID, names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::filter(!is.na(Intensity) & Intensity > 0)

  # Map LOD values
  lod_vec <- setNames(lod_df$LOD_Matrix, lod_df$Mediator)
  df_long$LOD <- lod_vec[df_long$Lipid]

  # Log transformation handling (applied to intensity and reference thresholds)
  x_label <- "Intensity (Raw)"
  if(isTRUE(config$missingness$log_transform_viz)) {
    df_long$Intensity <- log2(df_long$Intensity)
    df_long$LOD <- log2(df_long$LOD)
    x_label <- "Log2(Intensity)"
  }

  lod_lines <- df_long %>%
    dplyr::distinct(Lipid, LOD) %>%
    dplyr::filter(!is.na(LOD))

  # Histograms plot (including n=1 metrics natively)
  p <- ggplot2::ggplot(df_long) +
    ggplot2::geom_histogram(
      aes(x = Intensity),
      bins = config$missingness$hist_bins,
      fill = fill_color,
      color = "black",
      alpha = 0.7
    ) +
    ggplot2::geom_vline(
      data = lod_lines,
      aes(xintercept = LOD),
      color = "firebrick",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 7, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = x_label, y = "Count") +
    ggplot2::facet_wrap(~Lipid, scales = "free", ncol = 5)

  return(p)
}


# Worker: Plot Phenotype Completeness
#' Renders stacked completeness bar charts normalized by overall cohort size.
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
          # Prevent paste0() from coercing real NA values into the string "Quantile NA"
          !!binned_col_name := dplyr::if_else(
            is.na(!!rlang::sym(binned_col_name)),
            NA_character_,
            paste0("Quantile ", !!rlang::sym(binned_col_name))
          )
        )
      plot_col <- binned_col_name
    }

    # 3. Calculate and Plot
    # Filter out actual NA clinical entries
    df_valid <- df_join %>%
      dplyr::filter(!is.na(!!rlang::sym(plot_col)) & !!rlang::sym(plot_col) != "Quantile NA")

    # Track true cohort denominator to keep stacked heights <= 100%
    total_subjects <- length(unique(df_valid$Subject_ID))

    summary_df <- df_valid %>%
      dplyr::group_by(Lipid, !!rlang::sym(plot_col)) %>%
      dplyr::summarise(
        measured_count = sum(is_missing == 0),
        # Calculate contribution to cohort total
        comp_contrib = measured_count / total_subjects,
        .groups = "drop"
      )

    p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Lipid, y = comp_contrib, fill = as.factor(!!rlang::sym(plot_col)))) +
      ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top",
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = paste("Completeness Partitioned by", v_name),
        y = "Completeness (%)",
        x = "Lipid Species"
      )

    # 4. Dynamic Export
    file_name <- file.path(out_dir, paste0("pheno_completeness_", v_name, ".png"))
    ggplot2::ggsave(
      filename = file_name,
      plot = p,
      width = config$export$width,
      height = config$export$height,
      dpi = config$export$dpi,
      bg = "white"
    )

    plots[[v_name]] <- p
  }
  return(plots)
}


# PART 3: WORKERS (Statistical Engines - Stubs)----------------------------------

worker_run_fishers_exact <- function(...) {
  # TODO: Implement 2x2 contingency table math for binarized pipeline
}

worker_run_firth_logistic <- function(...) {
  # TODO: Implement logistf::logistf() for penalized logistic modeling
}

worker_run_tobit_model <- function(...) {
  # TODO: Implement survival::survreg() or AER::tobit() for left-censoring modeling
}

worker_run_hurdle_model <- function(...) {
  # TODO: Implement pscl::hurdle() modeling
}

worker_aggregate_lipid_classes <- function(...) {
  # TODO: Implement aggregation mapping to reconstruct functional pools
}


# SCRIPT: 00_lipid_engines.R
# PURPOSE: Master statistical engines and visualization wrappers for sparse lipidomics.
# ARCHITECTURE:
#   - Parts 1 & 2: Differential Expression (Binarized Hurdle 1 - Fisher/Firth)
#   - Parts 3 & 4: Raw Abundances (Continuous Hurdle 2 - Mann-Whitney/Spearman)
#   - Parts 5 & 6: Lipid Set Enrichment Analysis (LSEA Binary Burden Math)

library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
# library(logistf) # Required for Firth's Penalized Regression

# PART 4: CORE STATISTICAL WORKERS (DEA)-----------------------------------

#' Worker: Fisher's Exact Test (Categorical)
worker_fisher_exact <- function(v_valid, g_valid, target_name, ref_name) {
  # Build 2x2 table: rows = Group, cols = Detection (0, 1)
  tab <- table(factor(g_valid, levels = c(ref_name, target_name)),
               factor(v_valid, levels = c(0, 1)))

  if (all(dim(tab) == c(2, 2))) {
    res <- suppressWarnings(fisher.test(tab))
    return(list(p_val = res$p.value, estimate = res$estimate[[1]])) # Odds Ratio
  } else {
    return(list(p_val = NA_real_, estimate = NA_real_))
  }
}

#' Worker: Firth's Penalized Logistic Regression (Continuous)
worker_firth_regression <- function(v_valid, clin_vec_cont) {
  df <- data.frame(Detected = v_valid, Predictor = clin_vec_cont)
  # Uses logistf to handle complete separation in sparse omics
  fit <- try(logistf::logistf(Detected ~ Predictor, data = df), silent = TRUE)

  if (inherits(fit, "try-error")) return(list(p_val = NA_real_, estimate = NA_real_))

  p_val <- fit$prob["Predictor"]
  estimate <- exp(coef(fit)["Predictor"]) # Odds Ratio per 1-unit increase
  return(list(p_val = p_val, estimate = estimate))
}

# PART 5: DIFFERENTIAL EXPRESSION WRAPPER-----------------------------------------------------------

#' Wrapper: Sparse Differential Expression Analysis
wrap_sparse_dea <- function(mat, clin_df, config, engine, job_name, out_dir, assets_dir, project_colors_func) {

  # 1. Filter to requested contrast groups (if categorical)
  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_df <- clin_df %>% dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups))
  }

  common_samples <- intersect(colnames(mat), clin_df$Subject_ID)
  if (length(common_samples) < 5) return(list(status = "error", error_msg = "Insufficient subjects."))

  mat_sub <- mat[, common_samples, drop = FALSE]
  clin_sub <- clin_df %>% filter(Subject_ID %in% common_samples)
  clin_sub <- clin_sub[match(common_samples, clin_sub$Subject_ID), ]

  # Setup formatting
  target_name <- if(!is.null(config$target_groups)) paste(config$target_groups, collapse = "+") else "High"
  ref_name <- if(!is.null(config$ref_groups)) paste(config$ref_groups, collapse = "+") else "Low"

  if (engine == "fisher") {
    clin_sub$Plot_Group <- dplyr::case_when(
      clin_sub[[config$test_var]] %in% config$target_groups ~ target_name,
      clin_sub[[config$test_var]] %in% config$ref_groups ~ ref_name,
      TRUE ~ NA_character_
    )
    clin_vec <- clin_sub$Plot_Group
  } else {
    clin_vec <- clin_sub[[config$test_var]] # Continuous
  }

  # 2. Iterative Testing & Zero-Variance Exclusion
  res_list <- list()
  failed_list <- list()

  for (lipid in rownames(mat_sub)) {
    v_valid <- mat_sub[lipid, ]

    # Exclude if no variance (all 0s or all 1s in this subset)
    if (length(unique(v_valid)) <= 1) {
      failed_list[[lipid]] <- data.frame(Lipid = lipid, Method = "Zero Variance", stringsAsFactors = FALSE)
      next
    }

    # Route to appropriate engine
    if (engine == "fisher") {
      stat_res <- worker_fisher_exact(v_valid, clin_vec, target_name, ref_name)
    } else {
      stat_res <- worker_firth_regression(v_valid, clin_vec)
    }

    res_list[[lipid]] <- data.frame(
      Lipid = lipid,
      OddsRatio = stat_res$estimate,
      Log2OR = log2(stat_res$estimate + 1e-6), # Safe log2
      P_Value = stat_res$p_val,
      stringsAsFactors = FALSE
    )
  }

  # Combine results
  if (length(res_list) == 0) return(list(status = "error", error_msg = "All lipids failed variance checks."))

  full_res <- do.call(rbind, res_list)
  failed_res <- if(length(failed_list) > 0) do.call(rbind, failed_list) else data.frame()

  # 3. Benjamini-Hochberg FDR
  full_res$FDR <- p.adjust(full_res$P_Value, method = "BH")
  full_res <- full_res %>% arrange(FDR, P_Value)

  sig_res <- full_res %>% filter(FDR < config$dea_p_cutoff & abs(Log2OR) > config$dea_fc_cutoff)

  # 4. Generate Volcano Plot
  volcano_plot <- NULL
  if (nrow(full_res) > 0) {
    plot_df <- full_res %>%
      mutate(
        Sig = case_when(
          FDR < config$dea_p_cutoff & Log2OR > config$dea_fc_cutoff ~ "Up",
          FDR < config$dea_p_cutoff & Log2OR < -config$dea_fc_cutoff ~ "Down",
          TRUE ~ "NS"
        ),
        NegLog10P = -log10(P_Value)
      )

    volcano_plot <- ggplot(plot_df, aes(x = Log2OR, y = NegLog10P, color = Sig)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "grey80")) +
      geom_hline(yintercept = -log10(max(full_res$P_Value[full_res$FDR < config$dea_p_cutoff], na.rm=T)), linetype = "dashed") +
      geom_vline(xintercept = c(-config$dea_fc_cutoff, config$dea_fc_cutoff), linetype = "dashed") +
      theme_minimal() +
      labs(title = sprintf("Volcano: %s", config$title), x = "Log2 Odds Ratio", y = "-Log10 P-Value")

    ggsave(file.path(out_dir, paste0(job_name, "_volcano.png")), volcano_plot, width = 8, height = 6, bg = "white")
    if (!is.null(assets_dir)) ggsave(file.path(assets_dir, paste0(job_name, "_volcano.png")), volcano_plot, width = 8, height = 6, bg = "white")
  }

  # 5. Generate Heatmap (Top Hits)
  hm_plot <- NULL
  if (nrow(sig_res) > 2) {
    top_lipids <- head(sig_res$Lipid, 25)
    hm_mat <- mat_sub[top_lipids, ]

    col_fun <- colorRamp2(c(0, 1), c("white", "darkblue"))

    # Ensure safe fallback colors if passing categorical
    if (engine == "fisher") {
      base_levels <- c(config$target_groups[1], config$ref_groups[1])
      mapped_cols <- project_colors_func(base_levels)
      pheno_cols <- setNames(mapped_cols, c(target_name, ref_name))
    } else {
      pheno_cols <- circlize::colorRamp2(c(min(clin_vec), max(clin_vec)), c("white", "darkred"))
    }

    ha <- HeatmapAnnotation(Group = clin_vec, col = list(Group = pheno_cols))

    hm_plot <- Heatmap(hm_mat, name = "Detected", top_annotation = ha, col = col_fun,
                       show_row_names = TRUE, show_column_names = FALSE,
                       cluster_columns = TRUE, cluster_rows = TRUE,
                       column_title = config$title)

    png(file.path(out_dir, paste0(job_name, "_heatmap.png")), width = 8, height = 8, units = "in", res = 300)
    draw(hm_plot)
    invisible(dev.off())
    if (!is.null(assets_dir)) {
      png(file.path(assets_dir, paste0(job_name, "_heatmap.png")), width = 8, height = 8, units = "in", res = 300)
      draw(hm_plot)
      invisible(dev.off())
    }
  }

  return(list(
    status = "success",
    tables = list(full_results = full_res, significant_features = sig_res, failed_features = failed_res),
    plots = list(volcano = volcano_plot, heatmap = hm_plot)
  ))
}

# PART 6: RAW ABUNDANCES (Hurdle 2)----------------------------------------

#' Wrapper: Raw Abundance Plotting
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

  # 2. Statistical Testing (Conditional Abundance / Hurdle 2)
  stats_list <- lapply(rownames(mat_sub), function(lipid) {
    val_vec <- mat_sub[lipid, ]

    # Isolate strictly detected data
    valid_idx <- !is.na(val_vec) & val_vec > 0 & !is.na(clin_vec_plot)
    v_valid <- val_vec[valid_idx]
    g_valid <- clin_vec_plot[valid_idx]

    if (!is_cont) {
      # Mann-Whitney U for Categorical Contrasts
      v_target <- v_valid[g_valid == target_name]
      v_ref <- v_valid[g_valid == ref_name]
      if (length(v_target) >= 3 && length(v_ref) >= 3) {
        p_val <- suppressWarnings(wilcox.test(v_target, v_ref, exact = FALSE)$p.value)
      } else { p_val <- NA_real_ }
    } else {
      # Spearman Correlation for Continuous Contrasts
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
    grid_height <- max(4, ceiling(nrow(mat_sub) / 4) * 3.5)
    file_name <- paste0(job_name, "_raw_intensities.png")
    ggsave(file.path(out_dir, file_name), plot = p_raw, width = 16, height = grid_height, dpi = 300, bg = "white")
    if (!is.null(assets_dir)) ggsave(file.path(assets_dir, file_name), plot = p_raw, width = 16, height = grid_height, dpi = 300, bg = "white")
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

    # FIX: Pass both groups simultaneously to properly coordinate the fallback palette
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
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey30"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 10),
        panel.grid.minor = element_blank()
      ) +
      labs(title = lipid, subtitle = sub_str, y = "Log10(Intensity)")

    plot_list[[lipid]] <- p
  }

  if (length(plot_list) == 0) return(NULL)

  # 4 Column Patchwork Layout
  grid_plot <- patchwork::wrap_plots(plot_list, ncol = 4) +
    patchwork::plot_annotation(
      title = sprintf("Raw Conditional Abundances: %s", config$title),
      subtitle = "Values strictly log10-transformed after excluding non-detected features.",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  return(grid_plot)
}

# PART 7: LIPID SET ENRICHMENT ANALYSIS (LSEA) --------------------------------------

#' Wrapper: Lipid Set Enrichment Analysis (Binary Burden)
wrap_lsea <- function(burden_mat, clin_df, config, job_name, out_dir, assets_dir, project_colors_func) {

  if (!is.null(config$target_groups) && !is.null(config$ref_groups)) {
    clin_df <- clin_df %>% dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups))
  }

  common_samples <- intersect(colnames(burden_mat), clin_df$Subject_ID)
  if (length(common_samples) < 5) return(list(status = "error", error_msg = "Insufficient subjects for LSEA."))

  mat_sub <- burden_mat[, common_samples, drop = FALSE]
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

  # 2. Statistical Testing (Binary Burden)
  stats_list <- lapply(rownames(mat_sub), function(family) {
    v_valid <- mat_sub[family, ]
    g_valid <- clin_vec_plot

    # Skip families with zero variance
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

  write.csv(stats_df, file.path(out_dir, paste0(job_name, "_LSEA_results.csv")), row.names = FALSE)

  # 3. Call Plotting Worker
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

  # 4. Save dynamically scaled image
  if (!is.null(p_lsea)) {
    grid_height <- max(4, ceiling(nrow(mat_sub) / 4) * 4)
    file_name <- paste0(job_name, "_LSEA_burden.png")
    ggsave(file.path(out_dir, file_name), plot = p_lsea, width = 16, height = grid_height, dpi = 300, bg = "white")
    if (!is.null(assets_dir)) ggsave(file.path(assets_dir, file_name), plot = p_lsea, width = 16, height = grid_height, dpi = 300, bg = "white")
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

    # Safe Color Mapping
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
      # FIX 1: Set height = 0 to prevent vertical smearing into negative values
      geom_jitter(width = 0.2, height = 0, alpha = 0.8, size = 2) +
      scale_fill_manual(values = p_cols) +
      scale_color_manual(values = p_cols) +
      # FIX 2: Force Y-axis to strictly draw every single integer (step of 1)
      scale_y_continuous(breaks = function(limits) seq(0, ceiling(limits[2]), by = 1)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey30"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 10),
        panel.grid.minor = element_blank()
      ) +
      labs(title = family, subtitle = sub_str, y = "Detected Lipids (Count)")

    plot_list[[family]] <- p
  }

  if (length(plot_list) == 0) return(NULL)

  # 4 Column Patchwork Layout
  grid_plot <- patchwork::wrap_plots(plot_list, ncol = 4) +
    patchwork::plot_annotation(
      title = sprintf("LSEA Binary Burden: %s", config$title),
      subtitle = "Boxplots representing the total count of detected lipids per biological family.",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  return(grid_plot)
}

