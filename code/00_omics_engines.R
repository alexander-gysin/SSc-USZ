# ==============================================================================
# OMICS ENGINES: High-Dimensional Analysis Functions
# Project: Systemic Sclerosis Disease Activity Multi-Omics
# ==============================================================================

library(limma)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(ggrepel)
library(patchwork)
library(ggridges)

# ------------------------------------------------------------------------------
# 1. PCA ENGINE
# ------------------------------------------------------------------------------
run_pca_engine <- function(omics_mat, clin_df, pca_conf, project_colors) {

  # Step 1: Subset Data if requested
  if (!is.null(pca_conf$subset_to_contrast)) {
    target_col <- pca_conf$subset_to_contrast$group_col
    target_lvls <- pca_conf$subset_to_contrast$target_groups
    clin_df <- clin_df %>% filter(!!sym(target_col) %in% target_lvls)
  }

  valid_barcodes <- intersect(colnames(omics_mat), clin_df$Subject_ID)
  clin_df <- clin_df %>% filter(Subject_ID %in% valid_barcodes)
  mat_sub <- omics_mat[, clin_df$Subject_ID, drop = FALSE]

  # Step 2: Run PCA Math
  pca_res <- FactoMineR::PCA(t(mat_sub), scale.unit = TRUE, graph = FALSE)

  generated_plots <- list()

  # Step 3: Base Diagnostic Plots (NEW: Comprehensive Titles)
  generated_plots[["scree"]] <- fviz_eig(pca_res, addlabels = TRUE) +
    theme_minimal() + labs(title = sprintf("%s: Variance Explained (Scree)", pca_conf$title))

  generated_plots[["loadings"]] <- fviz_pca_var(pca_res, select.var = list(contrib = 30),
                                                col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
    theme_minimal() + labs(title = sprintf("%s: PCA Loadings (Top 30)", pca_conf$title))

  generated_plots[["pc1_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 1, top = 10, fill = "grey30", color = "grey30") +
    theme_minimal() + coord_flip() + labs(title = sprintf("%s: Top PC1 Drivers", pca_conf$title))

  generated_plots[["pc2_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 2, top = 10, fill = "grey30", color = "grey30") +
    theme_minimal() + coord_flip() + labs(title = sprintf("%s: Top PC2 Drivers", pca_conf$title))

  # Step 4: Dynamic Overlays based on Configuration
  for (i in seq_along(pca_conf$color_vars)) {
    c_var <- pca_conf$color_vars[i]

    if (!(c_var %in% colnames(clin_df))) {
      warning(sprintf("PCA Engine Warning: Column '%s' not found in clinical data. Skipping this overlay.", c_var))
      next
    }

    plot_idx <- paste0("overlay_", c_var)
    dynamic_title <- sprintf("%s (Colored by %s)", pca_conf$title, c_var) # NEW: Comprehensive title

    if (is.numeric(clin_df[[c_var]])) {
      p <- fviz_pca_ind(pca_res, label = "none", col.ind = clin_df[[c_var]], legend.title = c_var) +
        scale_color_viridis_c(option = "viridis") +
        theme_minimal() + theme(legend.position = "bottom") + labs(title = dynamic_title)
    } else {
      clin_df[[c_var]] <- as.character(clin_df[[c_var]])
      clin_df[[c_var]][is.na(clin_df[[c_var]])] <- "Missing"
      clin_df[[c_var]] <- as.factor(clin_df[[c_var]])

      dynamic_palette <- project_colors(levels(clin_df[[c_var]]))

      p <- fviz_pca_ind(pca_res, label = "none", habillage = clin_df[[c_var]], palette = dynamic_palette, addEllipses = TRUE) +
        theme_minimal() + theme(legend.position = "bottom") + labs(title = dynamic_title)
    }
    generated_plots[[plot_idx]] <- p
  }

  return(generated_plots)
}


# ------------------------------------------------------------------------------
# 2. DIFFERENTIAL EXPRESSION ENGINE (Limma)
# ------------------------------------------------------------------------------
run_limma_engine <- function(omics_mat, clin_df, dea_conf, omics_config, project_colors) {

  clin_sub <- clin_df %>% drop_na(all_of(c(dea_conf$group_col, dea_conf$covariates)))

  if (isTRUE(dea_conf$is_continuous)) {
    split_val <- if (dea_conf$split_method == "median") median(clin_sub[[dea_conf$group_col]], na.rm=TRUE) else mean(clin_sub[[dea_conf$group_col]], na.rm=TRUE)
    clin_sub <- clin_sub %>% mutate(limma_group = case_when(
      .data[[dea_conf$group_col]] > split_val ~ "Target",
      .data[[dea_conf$group_col]] <= split_val ~ "Reference",
      TRUE ~ NA_character_
    ))
  } else {
    clin_sub <- clin_sub %>% mutate(limma_group = case_when(
      .data[[dea_conf$group_col]] %in% dea_conf$target_groups ~ "Target",
      .data[[dea_conf$group_col]] %in% dea_conf$ref_groups ~ "Reference",
      TRUE ~ NA_character_
    ))
  }

  clin_sub <- clin_sub %>% filter(!is.na(limma_group)) %>% mutate(limma_group = factor(limma_group, levels = c("Reference", "Target")))
  if(length(unique(clin_sub$limma_group)) < 2) return(NULL)

  mat_sub <- omics_mat[, clin_sub$Subject_ID, drop = FALSE]

  formula_str <- if(length(dea_conf$covariates) > 0) paste("~ 0 + limma_group +", paste(dea_conf$covariates, collapse = " + ")) else "~ 0 + limma_group"

  design <- model.matrix(as.formula(formula_str), data = clin_sub)
  colnames(design)[1:2] <- c("Reference", "Target")

  fit <- lmFit(mat_sub, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Target - Reference, levels = design)))

  dea_res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P") %>%
    rownames_to_column("Feature") %>%
    mutate(
      Significance = case_when(
        adj.P.Val < omics_config$dea_p_cutoff & logFC > omics_config$dea_fc_cutoff  ~ "Up-regulated",
        adj.P.Val < omics_config$dea_p_cutoff & logFC < -omics_config$dea_fc_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      label_score = abs(logFC) * -log10(adj.P.Val) # NEW: Use adjusted P-value for labeling priority
    )

  top_labels <- bind_rows(
    dea_res %>% filter(Significance == "Up-regulated") %>% slice_max(order_by = label_score, n = 10),
    dea_res %>% filter(Significance == "Down-regulated") %>% slice_max(order_by = label_score, n = 10)
  )

  target_color <- if (!is.null(dea_conf$target_groups)) project_colors(dea_conf$target_groups[1]) else "#D55E00"
  ref_color <- if (!is.null(dea_conf$ref_groups)) project_colors(dea_conf$ref_groups[1]) else "#0072B2"

  # NEW: Plot using adj.P.Val and inject the FDR cutoff line
  p_volc <- ggplot(dea_res, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down-regulated" = ref_color, "Not Significant" = "grey85", "Up-regulated" = target_color)) +
    geom_text_repel(data = top_labels, aes(label = Feature), size = 3, color = "black", box.padding = 0.5, max.overlaps = Inf, show.legend = FALSE) +
    geom_vline(xintercept = c(-omics_config$dea_fc_cutoff, omics_config$dea_fc_cutoff), linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = -log10(omics_config$dea_p_cutoff), linetype = "dashed", color = "darkred") + # <--- FDR CUTOFF LINE
    theme_minimal() + theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold")) +
    labs(title = paste("Volcano:", dea_conf$title), x = "Log2 Fold Change", y = "-log10(FDR Adjusted P-value)")

  return(list(data = dea_res, volcano = p_volc, clin_used = clin_sub))
}

# ------------------------------------------------------------------------------
# 3. GSEA ENGINE
# ------------------------------------------------------------------------------
run_gsea_engine <- function(dea_res, go_db, omics_config, title) {

  gsea_input <- dea_res %>% filter(!is.na(Feature), !is.na(t)) %>%
    group_by(Feature) %>% slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(desc(t))

  ranked_vec <- setNames(gsea_input$t, gsea_input$Feature)[is.finite(gsea_input$t)]
  set.seed(42)

  gsea_res <- GSEA(geneList = ranked_vec, TERM2GENE = go_db, pvalueCutoff = 1,
                   minGSSize = omics_config$gsea_min_size, maxGSSize = 500, eps = 1e-10, nPermSimple = 10000, verbose = FALSE)

  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(NULL)

  gsea_res@result$Description <- gsub("GOBP_", "", gsea_res@result$Description) %>% gsub("_", " ", .)

  # NEW: Safe Filtering and Exploratory Fallback
  gsea_filtered <- as.data.frame(gsea_res) %>% filter(p.adjust < omics_config$gsea_p_cutoff) %>%
    mutate(Status = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>% arrange(desc(abs(NES)))

  warning_tag <- ""

  if(nrow(gsea_filtered) == 0) {
    warning_tag <- "\n(Exploratory: No paths passed strict FDR cutoff)"
    # Fallback to top 5 up and down by raw Enrichment Score
    gsea_filtered <- as.data.frame(gsea_res) %>%
      mutate(Status = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>%
      arrange(desc(abs(NES))) %>%
      group_by(Status) %>% slice_head(n = 5) %>% ungroup()
  }

  plot_df <- gsea_filtered %>% group_by(Status) %>% slice_max(order_by = abs(NES), n = 10) %>% ungroup() %>% mutate(Desc = str_wrap(Description, width = 40))

  p_dot <- ggplot(plot_df, aes(x = NES, y = reorder(Desc, NES))) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    scale_color_gradient(low = "firebrick3", high = "navy") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() + labs(title = paste0("GSEA: ", title, warning_tag), x = "Normalized Enrichment Score", y = NULL)

  sim_obj <- pairwise_termsim(gsea_res)

  # Prevent Ridgeplot from crashing if only a few exploratory pathways exist
  num_paths <- min(15, nrow(gsea_filtered))
  if(num_paths < 2) return(list(data = gsea_filtered, dotplot = p_dot, ridgeplot = NULL))

  p_ridge <- ridgeplot(sim_obj, showCategory = num_paths, core_enrichment = TRUE) +
    theme_minimal() + geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(title = paste0("Pathway Distribution: ", title, warning_tag), x = "t-statistic", y = NULL)

  return(list(data = gsea_filtered, dotplot = p_dot, ridgeplot = p_ridge))
}

# ------------------------------------------------------------------------------
# 4. CROSS-CONTRAST ENGINE (Trajectory Shifts)
# ------------------------------------------------------------------------------
run_cross_contrast_engine <- function(dea_x, dea_y, conf, omics_config, x_label, y_label, color_x, color_y, color_shared) {

  temporal_shifts <- full_join(
    dea_x %>% select(Feature, logFC_X = logFC, FDR_X = adj.P.Val),
    dea_y %>% select(Feature, logFC_Y = logFC, FDR_Y = adj.P.Val),
    by = "Feature"
  ) %>%
    mutate(
      Sig_X = FDR_X < omics_config$dea_p_cutoff & abs(logFC_X) > omics_config$dea_fc_cutoff,
      Sig_Y = FDR_Y < omics_config$dea_p_cutoff & abs(logFC_Y) > omics_config$dea_fc_cutoff,
      Signature = case_when(
        Sig_X & Sig_Y  ~ "Shared Drivers",
        Sig_X & !Sig_Y ~ "X-Axis Specific",
        !Sig_X & Sig_Y ~ "Y-Axis Specific",
        TRUE ~ "Not Significant"
      ),

      # 1. Pi-Score for ranking Shared Drivers globally
      label_score = case_when(
        Signature == "Shared Drivers"  ~ (abs(logFC_X) * -log10(FDR_X)) + (abs(logFC_Y) * -log10(FDR_Y)),
        TRUE ~ 0
      ),

      # 2. Discordance Score for ranking Specific Drivers (Distance from y=x)
      discordance_score = abs(logFC_X - logFC_Y)
    )

  # ONLY keep significant points. Background points are dropped completely.
  sig_points <- temporal_shifts %>% filter(Signature != "Not Significant")

  if(nrow(sig_points) == 0) return(NULL)

  # LABELING LOGIC
  top_n_spec <- if(!is.null(conf$top_n_specific)) conf$top_n_specific else 5
  top_n_shar <- if(!is.null(conf$top_n_shared)) conf$top_n_shared else 10

  # Grab Top Shared Drivers (Ranked globally by Combined Pi-Score)
  labels_shared <- sig_points %>%
    filter(Signature == "Shared Drivers") %>%
    slice_max(order_by = label_score, n = top_n_shar, with_ties = FALSE)

  # Grab Top Specific Drivers (Ranked by distance from the 45-degree line)
  labels_specific <- sig_points %>%
    filter(Signature %in% c("X-Axis Specific", "Y-Axis Specific")) %>%
    group_by(Signature) %>%
    slice_max(order_by = discordance_score, n = top_n_spec, with_ties = FALSE) %>%
    ungroup()

  top_labels <- bind_rows(labels_shared, labels_specific)

  # NEW: Obstacle-Aware Labeling.
  # Inject the labels back into the main dataset so ggrepel sees ALL points as obstacles
  sig_points <- sig_points %>%
    mutate(repel_label = ifelse(Feature %in% top_labels$Feature, Feature, ""))

  method_caption <- sprintf(
    "Labeling: Top %d Shared Drivers (Combined Pi-Score) + Top %d Specific Drivers furthest from 45° line.",
    top_n_shar, top_n_spec
  )

  p_shift <- ggplot() +
    # Layer 1: Gridlines
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = 0, color = "grey40") + geom_vline(xintercept = 0, color = "grey40") +

    # Layer 2: Bold Foreground Significant points ONLY
    geom_point(data = sig_points, aes(x = logFC_X, y = logFC_Y, color = Signature), alpha = 0.8, size = 2.5) +

    # Layer 3: Obstacle-Aware Text Labels
    geom_text_repel(
      data = sig_points,
      aes(x = logFC_X, y = logFC_Y, label = repel_label),
      color = "black",
      size = 3.5,
      box.padding = 0.8,       # Increased bumper around the text
      point.padding = 0.5,     # Increased bumper around the data point
      min.segment.length = 0,  # Always draw the connecting line
      max.overlaps = Inf,
      show.legend = FALSE
    ) +

    scale_color_manual(values = setNames(
      c(color_shared, color_x, color_y),
      c("Shared Drivers", "X-Axis Specific", "Y-Axis Specific")
    )) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face="bold"),
      plot.caption = element_text(hjust = 0, face = "italic", color = "grey40", size = 9, margin = margin(t = 15))
    ) +
    labs(
      title = conf$title,
      x = paste(x_label, "(Log2 Fold Change)"),
      y = paste(y_label, "(Log2 Fold Change)"),
      caption = method_caption
    ) + coord_fixed(ratio = 1)

  return(p_shift)
}
