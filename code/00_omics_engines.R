# ==============================================================================
# OMICS ENGINES: High-Dimensional Analysis Functions
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
library(ComplexHeatmap)
library(UpSetR)
library(glmnet)
library(circlize)
library(ggVennDiagram)

# ------------------------------------------------------------------------------
# 1. PCA PIPELINE WRAPPER
# ------------------------------------------------------------------------------
wrap_pca_pipeline <- function(omics_mat, clin_df, pca_conf, pca_id, base_output_dir, project_colors_func) {

  # 1. Setup Directories
  pca_base_dir <- file.path(base_output_dir, "PCA")
  pca_sub_dir <- file.path(pca_base_dir, pca_id)
  if (!dir.exists(pca_sub_dir)) dir.create(pca_sub_dir, recursive = TRUE)

  # 2. Data Filtering & Subset
  if (!is.null(pca_conf$subset_to_contrast)) {
    target_col <- pca_conf$subset_to_contrast$group_col
    target_lvls <- pca_conf$subset_to_contrast$target_groups
    clin_df <- clin_df %>% filter(!!sym(target_col) %in% target_lvls)
  }

  valid_barcodes <- intersect(colnames(omics_mat), clin_df$Subject_ID)
  clin_df <- clin_df %>% filter(Subject_ID %in% valid_barcodes)
  mat_sub <- omics_mat[, clin_df$Subject_ID, drop = FALSE]

  # 3. PCA Calculation
  pca_res <- FactoMineR::PCA(t(mat_sub), scale.unit = TRUE, graph = FALSE)
  plots <- list()

  # 4. Generate Individual Plots
  plots[["scree"]] <- fviz_eig(pca_res, addlabels = TRUE) +
    theme_project_base() +
    labs(title = sprintf("%s: Variance Explained (Scree)", pca_conf$title))

  # ADJUSTMENT: Reset barheight and increased barwidth to fix horizontal label collision
  plots[["loadings"]] <- fviz_pca_var(pca_res, select.var = list(contrib = 30),
                                      col.var = "contrib",
                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                      repel = TRUE) +
    theme_project_base() +
    guides(colour = guide_colorbar(barwidth = unit(8, "cm"), barheight = unit(0.3, "cm"))) +
    labs(title = sprintf("%s: PCA Loadings (Top 30)", pca_conf$title))

  plots[["pc1_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 1, top = 10, fill = "grey30", color = "grey30") +
    theme_project_base() + coord_flip() + labs(title = "Top PC1 Drivers")

  plots[["pc2_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 2, top = 10, fill = "grey30", color = "grey30") +
    theme_project_base() + coord_flip() + labs(title = "Top PC2 Drivers")

  # 5. Variable Overlays (Individual and Merged)
  overlay_plots_for_grid <- list()
  for (c_var in pca_conf$color_vars) {
    if (!(c_var %in% colnames(clin_df))) next

    dynamic_title <- sprintf("%s (Colored by %s)", pca_conf$title, c_var)

    if (is.numeric(clin_df[[c_var]])) {
      p <- fviz_pca_ind(pca_res, label = "none", col.ind = clin_df[[c_var]], legend.title = c_var) +
        scale_color_viridis_c(option = "viridis") + theme_project_base() + labs(title = dynamic_title)
    } else {
      clin_df[[c_var]] <- as.factor(clin_df[[c_var]])
      p <- fviz_pca_ind(pca_res, label = "none", habillage = clin_df[[c_var]],
                        palette = project_colors_func(levels(clin_df[[c_var]])),
                        addEllipses = TRUE) + theme_project_base() + labs(title = dynamic_title)
    }

    # Return individual overlay plots with dynamic names
    plots[[paste0("PCA_Plot_", c_var)]] <- p
    overlay_plots_for_grid[[c_var]] <- p + labs(title = sprintf("Colored by %s", c_var))
  }

  # 6. Automated Saving (Individual Plots)
  for (p_name in names(plots)) {
    ggsave(filename = file.path(pca_sub_dir, paste0(pca_id, "_", p_name, ".png")),
           plot = plots[[p_name]], width = 7, height = 6, dpi = 300, bg = "white")
  }

  # 7. Stitch & Save Merged Drivers
  driver_grid <- (plots[["pc1_drivers"]] | plots[["pc2_drivers"]]) +
    plot_annotation(title = sprintf("%s: PCA Variance Drivers", pca_conf$title),
                    tag_levels = 'A', theme = theme(plot.title = element_text(size = 16, face = "bold")))

  ggsave(file.path(pca_sub_dir, paste0(pca_id, "_Merged_Drivers.png")), driver_grid, width = 12, height = 6, dpi = 300, bg="white")
  plots[["merged_drivers"]] <- driver_grid

  # 8. Stitch & Save Merged Overlays
  if (length(overlay_plots_for_grid) > 0) {
    ncol_val <- if(length(overlay_plots_for_grid) > 2) 2 else 1
    master_scatter_grid <- wrap_plots(overlay_plots_for_grid, ncol = ncol_val) +
      plot_annotation(title = pca_conf$title, subtitle = "Colored Scatter Overlays",
                      tag_levels = 'A', theme = theme(plot.title = element_text(size = 18, face = "bold")))

    ggsave(file.path(pca_sub_dir, paste0(pca_id, "_Merged_Scatter_Overlays.png")), master_scatter_grid,
           width = 14, height = (6 * ceiling(length(overlay_plots_for_grid)/2)), dpi = 300, bg="white")

    plots[["merged_overlays"]] <- master_scatter_grid
  }

  return(plots)
}

# ------------------------------------------------------------------------------
# 2. DEA & DISCOVERY PIPELINE WRAPPER (Updated with Ridgeplot Saving)
# ------------------------------------------------------------------------------
wrap_dea_pipeline <- function(omics_mat, clin_df, dea_conf, dea_id, base_output_dir,
                              go_db, omics_config, heatmap_config, project_colors_func) {

  results <- list(status = "success", plots = list(), tables = list(), data = list())

  # 1. Setup Directories
  dea_base_dir <- file.path(base_output_dir, "DEA")
  dea_sub_dir <- file.path(dea_base_dir, dea_id)
  if (!dir.exists(dea_sub_dir)) dir.create(dea_sub_dir, recursive = TRUE)

  # 2. Pre-flight: Diagnostic Filtering
  if("Sample_Age" %in% colnames(clin_df)) clin_df$Sample_Age <- as.numeric(clin_df$Sample_Age)

  req_cols <- c(dea_conf$group_col, dea_conf$covariates)
  available_cols <- intersect(req_cols, colnames(clin_df))

  clin_sub <- clin_df %>%
    filter(!!sym(dea_conf$group_col) %in% c(dea_conf$target_groups, dea_conf$ref_groups)) %>%
    drop_na(all_of(available_cols))

  n_target <- sum(clin_sub[[dea_conf$group_col]] %in% dea_conf$target_groups)
  n_ref <- sum(clin_sub[[dea_conf$group_col]] %in% dea_conf$ref_groups)
  results$diagnostics <- list(n_target = n_target, n_ref = n_ref)

  if(n_target < 2 || n_ref < 2) {
    results$status <- "failed"
    results$error_msg <- sprintf("Too few samples after NA filter (%d T vs %d R)", n_target, n_ref)
    return(results)
  }

  # 3. Execution Phase ---------------------------------------------------------

  # A. Limma DEA
  limma_res <- run_limma_engine(omics_mat, clin_df, dea_conf, project_colors_func)
  if (is.null(limma_res)) {
    results$status <- "failed"; results$error_msg <- "Limma engine failed"; return(results)
  }
  results$data$dea <- limma_res$data
  results$plots$volcano <- limma_res$volcano

  write_csv(limma_res$data, file.path(dea_sub_dir, paste0(dea_id, "_DEA_Full_Results.csv")))
  ggsave(file.path(dea_sub_dir, paste0(dea_id, "_Volcano.png")), limma_res$volcano, width = 8, height = 7, dpi=300)

  # B. Significant Table extraction
  sig_table <- limma_res$data %>%
    filter(Significance != "Not Significant") %>%
    select(Feature, logFC, adj.P.Val, Significance) %>%
    arrange(adj.P.Val)

  results$tables$significant_proteins <- sig_table
  write_csv(sig_table, file.path(dea_sub_dir, paste0(dea_id, "_Significant_Proteins.csv")))

  # C. Heatmap
  results$plots$heatmap <- run_heatmap_engine(omics_mat, clin_sub, limma_res$data,
                                              title = dea_conf$title,
                                              n_top = heatmap_config$top_n,
                                              split_by_group = heatmap_config$split_by_group)

  png(file.path(dea_sub_dir, paste0(dea_id, "_Heatmap_Top", heatmap_config$top_n, ".png")),
      width=10, height=12, units="in", res=300)
  draw(results$plots$heatmap); dev.off()

  # D. GSEA
  gsea_res <- run_gsea_engine(limma_res$data, go_db, omics_config, dea_conf$title)
  if (!is.null(gsea_res)) {
    results$data$gsea <- gsea_res$data
    results$plots$gsea_dotplot <- gsea_res$dotplot
    results$plots$gsea_ridge <- gsea_res$ridgeplot

    write_csv(gsea_res$data, file.path(dea_sub_dir, paste0(dea_id, "_GSEA_Results.csv")))
    ggsave(file.path(dea_sub_dir, paste0(dea_id, "_GSEA_Dotplot.png")), gsea_res$dotplot, width = 9, height = 7, dpi=300)
    if(!is.null(gsea_res$ridgeplot)) {
      ggsave(file.path(dea_sub_dir, paste0(dea_id, "_GSEA_Ridgeplot.png")), gsea_res$ridgeplot, width = 9, height = 8, dpi=300)
    }
  }

  # E. LASSO (RESTORED CSV EXPORT)
  lasso_res <- run_lasso_engine(omics_mat, clin_df, dea_conf$group_col,
                                c(dea_conf$target_groups, dea_conf$ref_groups), title = dea_conf$title)
  if (!is.null(lasso_res)) {
    results$data$lasso <- lasso_res$data
    results$plots$lasso_cv <- lasso_res$cv_plot
    results$plots$lasso_weights <- lasso_res$coeff_plot

    # RESTORED: Exporting coefficients to CSV
    write_csv(lasso_res$data, file.path(dea_sub_dir, paste0(dea_id, "_LASSO_Coefficients.csv")))

    png(file.path(dea_sub_dir, paste0(dea_id, "_LASSO_CV.png")), width=7, height=6, units="in", res=300)
    plot(lasso_res$cv_plot); title(main = sprintf("%s: LASSO Tuning", dea_conf$title), line = 2.5); dev.off()

    ggsave(file.path(dea_sub_dir, paste0(dea_id, "_LASSO_Weights.png")), lasso_res$coeff_plot, width=7, height=6, dpi=300)
  }

  return(results)
}

# ------------------------------------------------------------------------------
# 3. TRAJECTORY MAPPING PIPELINE WRAPPER
# ------------------------------------------------------------------------------
wrap_trajectory_pipeline <- function(cross_conf, cross_id, master_dea, master_gsea,
                                     contrast_configs, base_output_dir,
                                     project_colors_func, highlight_color, omics_config) {

  results <- list(plots = list(), status = "success")

  # 1. Fetch Source Data & Setup Labels
  dea_x <- master_dea[[cross_conf$x_axis_contrast]]
  dea_y <- master_dea[[cross_conf$y_axis_contrast]]
  gsea_x <- master_gsea[[cross_conf$x_axis_contrast]]
  gsea_y <- master_gsea[[cross_conf$y_axis_contrast]]

  if (is.null(dea_x) || is.null(dea_y)) {
    results$status <- "failed"; results$error_msg <- "Source DEA results missing"; return(results)
  }

  name_x <- contrast_configs[[cross_conf$x_axis_contrast]]$title
  name_y <- contrast_configs[[cross_conf$y_axis_contrast]]$title
  vs_label <- sprintf("%s VS %s", name_x, name_y)

  color_x <- project_colors_func(contrast_configs[[cross_conf$x_axis_contrast]]$target_groups[1])
  color_y <- project_colors_func(contrast_configs[[cross_conf$y_axis_contrast]]$target_groups[1])

  traj_sub_dir <- file.path(base_output_dir, "Trajectories", cross_id)
  if (!dir.exists(traj_sub_dir)) dir.create(traj_sub_dir, recursive = TRUE)

  # ----------------------------------------------------------------------------
  # A. PROTEIN TRAJECTORY
  # ----------------------------------------------------------------------------
  # Logic: Merge DEAs and rank by pi-score
  results$plots$protein_trajectory <- run_cross_contrast_engine(
    dea_x, dea_y, cross_conf, name_x, name_y, color_x, color_y, highlight_color,
    title = paste("Protein Trajectory:", vs_label)
  )
  if(!is.null(results$plots$protein_trajectory)) {
    ggsave(file.path(traj_sub_dir, paste0(cross_id, "_Protein_Trajectory.png")),
           results$plots$protein_trajectory, width=9, height=9, dpi=300, bg="white")
  }

  # ----------------------------------------------------------------------------
  # B. VENN DIAGRAM (Logic Internalized)
  # ----------------------------------------------------------------------------
  sig_x <- dea_x %>% filter(Significance != "Not Significant") %>% pull(Feature)
  sig_y <- dea_y %>% filter(Significance != "Not Significant") %>% pull(Feature)

  venn_data <- list(only_x = length(setdiff(sig_x, sig_y)),
                    only_y = length(setdiff(sig_y, sig_x)),
                    shared = length(intersect(sig_x, sig_y)))

  if(sum(unlist(venn_data)) > 0) {
    theta <- seq(0, 2 * pi, length.out = 100)
    p_venn <- ggplot() +
      geom_polygon(data = data.frame(x = -0.5 + cos(theta), y = sin(theta)), aes(x, y), fill = color_x, alpha = 0.3, color = color_x, linewidth = 1.2) +
      geom_polygon(data = data.frame(x = 0.5 + cos(theta), y = sin(theta)), aes(x, y), fill = color_y, alpha = 0.3, color = color_y, linewidth = 1.2) +
      annotate("text", x = -0.6, y = 1.3, label = name_x, size = 4.5, fontface = "bold", color = color_x) +
      annotate("text", x = 0.6, y = 1.3, label = name_y, size = 4.5, fontface = "bold", color = color_y) +
      annotate("text", x = -0.9, y = 0, label = venn_data$only_x, size = 6, fontface = "bold") +
      annotate("text", x = 0.9, y = 0, label = venn_data$only_y, size = 6, fontface = "bold") +
      annotate("text", x = 0, y = 0, label = venn_data$shared, size = 7, fontface = "bold") +
      theme_void() + coord_fixed(clip = "off", xlim = c(-1.8, 1.8), ylim = c(-1.5, 1.5)) +
      labs(title = paste("DEA Venn:", vs_label)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.4)))

    results$plots$venn <- p_venn
    ggsave(file.path(traj_sub_dir, paste0(cross_id, "_Venn_Intersection.png")), p_venn, width=8, height=6, dpi=300, bg="white")
  }

  # ----------------------------------------------------------------------------
  # C. PATHWAY TRAJECTORY (Upgraded Selection & Coloring Logic)
  # ----------------------------------------------------------------------------
  if (!is.null(gsea_x) && !is.null(gsea_y)) {

    # 1. Merge and Define Quadrants/Significance
    p_cutoff <- omics_config$gsea_p_cutoff
    rank_var <- if(!is.null(cross_conf$pathway_rank_var)) cross_conf$pathway_rank_var else "p.adjust"
    n_per_quad <- if(!is.null(cross_conf$top_n_pathways)) cross_conf$top_n_pathways else 10

    merged_gsea <- inner_join(
      gsea_x %>% select(Description, NES_X = NES, p_X = p.adjust),
      gsea_y %>% select(Description, NES_Y = NES, p_Y = p.adjust),
      by = "Description"
    ) %>% mutate(
      Sig_X = p_X < p_cutoff, Sig_Y = p_Y < p_cutoff,
      Signature = case_when(Sig_X & Sig_Y ~ "Shared Drivers", Sig_X & !Sig_Y ~ sprintf("%s Specific", name_x), !Sig_X & Sig_Y ~ sprintf("%s Specific", name_y), TRUE ~ "Not Significant"),
      Quadrant = case_when(NES_X > 0 & NES_Y > 0 ~ "Up-Up", NES_X < 0 & NES_Y < 0 ~ "Down-Down", NES_X > 0 & NES_Y < 0 ~ "Up-Down", NES_X < 0 & NES_Y > 0 ~ "Down-Up", TRUE ~ "None"),
      rank_score = (p_X + p_Y) # Lowest combined FDR
    )

    # 2. Select Top N per Quadrant
    top_paths <- merged_gsea %>%
      filter(Signature != "Not Significant") %>%
      group_by(Quadrant) %>%
      slice_min(order_by = !!sym(rank_var), n = n_per_quad, with_ties = FALSE) %>%
      ungroup()

    # 3. Generate Plot
    method_caption <- sprintf("Top %d pathways per quadrant selected by lowest %s (Significant in at least one contrast).", n_per_quad, rank_var)
    color_map <- setNames(c("black", color_x, color_y), c("Shared Drivers", sprintf("%s Specific", name_x), sprintf("%s Specific", name_y)))

    p_path <- ggplot(merged_gsea, aes(x = NES_X, y = NES_Y)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
      geom_vline(xintercept = 0, color = "grey40") + geom_hline(yintercept = 0, color = "grey40") +
      geom_point(aes(color = Signature), alpha = 0.5, size = 2) +
      geom_text_repel(data = top_paths, aes(label = Description), size = 3, color = "black", box.padding = 0.5, max.overlaps = Inf) +
      scale_color_manual(values = color_map) +
      theme_project_base() +
      labs(title = str_wrap(paste("Pathway Trajectory:", vs_label), 60),
           x = paste(name_x, "(NES)"), y = paste(name_y, "(NES)"),
           caption = method_caption)

    results$plots$pathway_trajectory <- p_path
    ggsave(file.path(traj_sub_dir, paste0(cross_id, "_Pathway_Trajectory.png")), p_path, width=11, height=9, dpi=300, bg="white")
  }

  return(results)
}
