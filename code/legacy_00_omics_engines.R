# OMICS ENGINES: High-Dimensional Analysis Functions ---------------------------

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
library(variancePartition)
library(BiocParallel)
library(future)
library(future.apply)
library(car)

# PART 1: THE MASTER PIPELINE WRAPPERS (The "Managers") ------------------------

# 1. PCA Pipeline Wrapper ------------------------------------------------------
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

# 2. DEA & Discovery Pipeline Wrapper ------------------------------------------
# Description:
#   Orchestrates the entire Differential Expression Analysis (DEA) workflow.
#   It cleans the data, runs the Limma engine for significance, generates Volcano
#   plots, extracts top features for a Heatmap, and performs GSEA pathway analysis.
# Args:
#   omics_mat: Cleaned omics matrix.
#   clin_df: Clinical metadata dataframe.
#   dea_conf: List with contrast details (groups, covariates, cutoffs).
#   dea_id: String identifier for this run.
#   base_output_dir: Path for saving outputs.
#   go_db: Gene Ontology/Pathway database object for GSEA.
#   omics_config: Global config settings (p-value cutoffs, min set size).
#   heatmap_config: Settings for heatmap rendering.
#   project_colors_func: Function for cohort color mapping.
# Returns:
#   A nested list containing 'status', 'plots', 'tables', and 'data'.
wrap_dea_pipeline <- function(omics_mat, clin_df, dea_conf, dea_id, base_output_dir,
                              go_db, omics_config, heatmap_config, project_colors_func,
                              assets_dir = NULL) { # <-- NEW: Defaults to NULL for legacy support

  results <- list(status = "success", plots = list(), tables = list(), data = list())

  # Setup Directories
  dea_base_dir <- file.path(base_output_dir, "DEA")
  dea_sub_dir <- file.path(dea_base_dir, dea_id)
  if (!dir.exists(dea_sub_dir)) dir.create(dea_sub_dir, recursive = TRUE)

  # Ensure assets directory exists if provided
  if (!is.null(assets_dir) && !dir.exists(assets_dir)) {
    dir.create(assets_dir, recursive = TRUE)
  }

  # Pre-flight: Diagnostic Filtering
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

  # Execution Phase
  # A. Limma DEA
  limma_res <- run_limma_engine(omics_mat, clin_df, dea_conf, project_colors_func)
  if (is.null(limma_res)) {
    results$status <- "failed"; results$error_msg <- "Limma engine failed"; return(results)
  }
  results$data$dea <- limma_res$data
  results$plots$volcano <- limma_res$volcano

  # Save Legacy
  write_csv(limma_res$data, file.path(dea_sub_dir, paste0(dea_id, "_DEA_Full_Results.csv")))
  ggsave(file.path(dea_sub_dir, paste0(dea_id, "_Volcano.png")), limma_res$volcano, width = 8, height = 7, dpi=300)

  # Dual-Save Asset (normalized name)
  if (!is.null(assets_dir)) {
    ggsave(file.path(assets_dir, paste0(dea_id, "_volcano.png")), limma_res$volcano, width = 8, height = 7, dpi=300, bg="white")
  }

  # B. Significant Table extraction
  sig_table <- limma_res$data %>%
    filter(Significance != "Not Significant") %>%
    select(Feature, logFC, adj.P.Val, Significance) %>%
    arrange(adj.P.Val)

  results$tables$significant_proteins <- sig_table
  write_csv(sig_table, file.path(dea_sub_dir, paste0(dea_id, "_Significant_Proteins.csv")))

  # C. Heatmap
  custom_map <- if(!is.null(dea_conf$custom_color_map)) dea_conf$custom_color_map else NULL

  results$plots$heatmap <- run_heatmap_engine(omics_mat, clin_sub, limma_res$data,
                                              title = dea_conf$title,
                                              n_top = heatmap_config$top_n,
                                              split_by_group = heatmap_config$split_by_group,
                                              split_col = dea_conf$group_col,
                                              custom_color_map = custom_map,
                                              annotation_vars = heatmap_config$annotation_vars)

  # Save Legacy
  png(file.path(dea_sub_dir, paste0(dea_id, "_Heatmap_Top", heatmap_config$top_n, ".png")),
      width=10, height=12, units="in", res=300)
  ComplexHeatmap::draw(results$plots$heatmap)
  invisible(dev.off())

  # Dual-Save Asset (normalized name)
  if (!is.null(assets_dir)) {
    png(file.path(assets_dir, paste0(dea_id, "_heatmap.png")), width=10, height=12, units="in", res=300)
    ComplexHeatmap::draw(results$plots$heatmap)
    invisible(dev.off())
  }

  # D. Multi-Database GSEA
  results$plots$gsea <- list()
  results$data$gsea <- list()

  for (db_name in names(pathway_dbs)) {
    db_obj <- pathway_dbs[[db_name]]

    # Run the engine for this specific database
    gsea_res <- run_gsea_engine(limma_res$data, db_obj, omics_config, dea_conf$title, db_name)

    if (!is.null(gsea_res)) {
      results$data$gsea[[db_name]] <- gsea_res$data
      results$plots$gsea[[db_name]] <- gsea_res

      # Save Legacy CSV
      write_csv(gsea_res$data, file.path(dea_sub_dir, paste0(dea_id, "_", db_name, "_GSEA_Results.csv")))

      # Plot 1: Dotplot
      ggsave(file.path(dea_sub_dir, paste0(dea_id, "_", db_name, "_GSEA_Dotplot.png")), gsea_res$dotplot, width = 9, height = 7, dpi=300, bg="white")
      if (!is.null(assets_dir)) {
        ggsave(file.path(assets_dir, paste0(dea_id, "_", db_name, "_dotplot.png")), gsea_res$dotplot, width = 9, height = 7, dpi=300, bg="white")
      }

      # Plot 2: Ridgeplot
      if(!is.null(gsea_res$ridgeplot)) {
        ggsave(file.path(dea_sub_dir, paste0(dea_id, "_", db_name, "_GSEA_Ridgeplot.png")), gsea_res$ridgeplot, width = 9, height = 8, dpi=300, bg="white")
        if (!is.null(assets_dir)) {
          ggsave(file.path(assets_dir, paste0(dea_id, "_", db_name, "_ridgeplot.png")), gsea_res$ridgeplot, width = 9, height = 8, dpi=300, bg="white")
        }
      }

      # Plot 3: Emap
      if(!is.null(gsea_res$emap)) {
        ggsave(file.path(dea_sub_dir, paste0(dea_id, "_", db_name, "_GSEA_Emap.png")), gsea_res$emap, width = 10, height = 10, dpi=300, bg="white")
        if (!is.null(assets_dir)) {
          ggsave(file.path(assets_dir, paste0(dea_id, "_", db_name, "_emap.png")), gsea_res$emap, width = 10, height = 10, dpi=300, bg="white")
        }
      }
    }
  }

  return(results)
}

# 3. Trajectory Mapping Pipeline Wrapper ---------------------------------------
# Description:
#   Takes the outputs of two different DEA contrasts (e.g., A vs C and B vs C)
#   and intersects them to find shared and specific biological drivers. Generates
#   4-quadrant scatter plots for proteins and pathways, plus a Venn diagram.
# Args:
#   cross_conf: Config specifying which two contrasts to compare.
#   cross_id: String identifier for this trajectory analysis.
#   master_dea: Named list of previously generated DEA results.
#   master_gsea: Named list of previously generated GSEA results.
#   contrast_configs: Metadata definitions of the original contrasts.
#   base_output_dir: Root directory for saving results.
#   project_colors_func: Color mapping function.
#   highlight_color: Color to use for "Shared" significant features.
#   omics_config: Global config settings.
# Returns:
#   A list of combined multi-contrast plots (trajectory, venn, pathways).
wrap_trajectory_pipeline <- function(cross_conf, cross_id, master_dea, master_gsea,
                                     contrast_configs, base_output_dir,
                                     project_colors_func, highlight_color, omics_config) {

  results <- list(plots = list(), status = "success")

  # Fetch Source Data & Setup Labels
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

  # A. PROTEIN TRAJECTORY
  results$plots$protein_trajectory <- run_cross_contrast_engine(
    dea_x, dea_y, cross_conf, name_x, name_y, color_x, color_y, highlight_color,
    title = paste("Protein Trajectory:", vs_label)
  )
  if(!is.null(results$plots$protein_trajectory)) {
    ggsave(file.path(traj_sub_dir, paste0(cross_id, "_Protein_Trajectory.png")),
           results$plots$protein_trajectory, width=9, height=9, dpi=300, bg="white")
  }

  # B. VENN DIAGRAM
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

  # C. PATHWAY TRAJECTORY
  if (!is.null(gsea_x) && !is.null(gsea_y)) {

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
      rank_score = (p_X + p_Y)
    )

    top_paths <- merged_gsea %>%
      filter(Signature != "Not Significant") %>%
      group_by(Quadrant) %>%
      slice_min(order_by = !!sym(rank_var), n = n_per_quad, with_ties = FALSE) %>%
      ungroup()

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

# 4. Predictive Modeling Pipeline Wrapper --------------------------------------
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

# 5a. Joint VPA Wrapper --------------------------------------------------------
wrap_vpa_pipeline <- function(omics_mat, clin_df, vpa_conf, vpa_id, base_output_dir, top_n_plot = 30) {

  results <- list(status = "success", plots = list(), data = list())

  vpa_sub_dir <- file.path(base_output_dir, "VPA_Joint", vpa_id)
  if (!dir.exists(vpa_sub_dir)) dir.create(vpa_sub_dir, recursive = TRUE)

  min_minority <- if(!is.null(vpa_conf$min_minority_pct)) vpa_conf$min_minority_pct else 0.10

  # 1. Variable Safety Checks
  existing_vars <- intersect(vpa_conf$variables, colnames(clin_df))
  if(length(existing_vars) == 0) return(list(status = "failed", error_msg = "No matching variables found."))

  # 2. Align Data, Purge NAs, and Scale Numerics
  clin_sub <- clin_df %>%
    select(Subject_ID, any_of(existing_vars)) %>%
    drop_na() %>%
    droplevels() %>%
    mutate(across(where(is.numeric), ~ as.numeric(scale(.))))

  valid_ids <- intersect(colnames(omics_mat), clin_sub$Subject_ID)
  clin_sub <- clin_sub %>% filter(Subject_ID %in% valid_ids)

  if (nrow(clin_sub) < 15) {
    return(list(status = "failed", error_msg = sprintf("NA-Wipeout: Only %d patients remained after requiring complete cases.", nrow(clin_sub))))
  }

  # 3. Build Formula & Protect Against Ghost Levels + Class Imbalance
  formula_terms <- sapply(existing_vars, function(v) {
    if ((is.factor(clin_sub[[v]]) || is.character(clin_sub[[v]])) && length(unique(clin_sub[[v]])) < 2) {
      return(NULL)
    }

    is_categorical <- is.factor(clin_sub[[v]]) || is.character(clin_sub[[v]]) || is.logical(clin_sub[[v]])
    if (is_categorical) {
      tbl <- prop.table(table(clin_sub[[v]]))
      if (length(tbl) < 2 || min(tbl) < min_minority) return(NULL)
    }

    num_levels <- length(unique(clin_sub[[v]]))

    if (is_categorical && num_levels >= 5) {
      paste0("(1|`", v, "`)")
    } else {
      paste0("`", v, "`")
    }
  })

  formula_terms <- formula_terms[!sapply(formula_terms, is.null)]
  if(length(formula_terms) == 0) return(list(status="failed", error_msg="All variables had zero variance or extreme imbalance."))

  vpa_formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))

  # 4. Execution (PARALLELIZED ACROSS 4 CORES)
  res_varPart <- tryCatch({
    fitExtractVarPartModel(omics_mat[, valid_ids], vpa_formula, clin_sub, BPPARAM = SnowParam(workers = 4))
  }, error = function(e) return(e$message))

  if (is.character(res_varPart)) {
    return(list(status = "failed", error_msg = res_varPart))
  }

  # 5. Process & Plot
  varPart_df <- as.data.frame(res_varPart) %>% rownames_to_column("Feature")
  summary_df <- data.frame(
    Covariate = colnames(res_varPart),
    Mean_Var = colMeans(res_varPart, na.rm = TRUE),
    SEM = apply(res_varPart, 2, sd, na.rm = TRUE) / sqrt(nrow(res_varPart))
  ) %>%
    mutate(CI_upper = pmin(1, Mean_Var + (1.96 * SEM))) %>%
    arrange(desc(Mean_Var))

  plot_summary <- head(summary_df, top_n_plot)

  p_bar <- ggplot(plot_summary, aes(x = reorder(Covariate, Mean_Var), y = Mean_Var, fill = Covariate)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, Mean_Var - (1.96 * SEM)), ymax = CI_upper), width = 0.2) +
    geom_text(aes(y = CI_upper, label = sprintf("%.1f%%", Mean_Var * 100)), hjust = -0.2, fontface = "bold") +
    coord_flip() + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +
    theme_project_base() + theme(legend.position = "none") +
    labs(title = vpa_conf$title, subtitle = "Average Variance Explained", y = "% Variance Explained (± 95% CI)", x = NULL)

  long_var <- varPart_df %>%
    select(Feature, all_of(plot_summary$Covariate)) %>%
    pivot_longer(-Feature, names_to = "Covariate", values_to = "Val") %>%
    mutate(Covariate = factor(Covariate, levels = rev(plot_summary$Covariate)))

  p_violin <- ggplot(long_var, aes(x = Covariate, y = Val, fill = Covariate)) +
    geom_violin(alpha = 0.6, color = "grey40", scale = "width") +
    geom_boxplot(width = 0.08, fill = "white", color = "black", outlier.shape = 19, outlier.size = 1, outlier.alpha = 0.4) +
    coord_flip() + scale_y_continuous(labels = scales::percent) +
    theme_project_base() + theme(legend.position = "none") +
    labs(title = vpa_conf$title, x = NULL, y = "% Variance Explained")

  # 6. Package Results
  results$plots <- list(bar = p_bar, violin = p_violin)
  results$data <- varPart_df
  readr::write_csv(varPart_df, file.path(vpa_sub_dir, paste0(vpa_id, "_Joint_VPA_Results.csv")))
  ggsave(file.path(vpa_sub_dir, paste0(vpa_id, "_Bar.png")), p_bar, width = 8, height = max(6, nrow(plot_summary)*0.3), bg="white")

  return(results)
}

# 5b. Univariate VPA Scan Wrapper ----------------------------------------------
wrap_vpa_univariate_pipeline <- function(omics_mat, clin_df, vpa_conf, vpa_id, base_output_dir) {

  # Setup Directories
  vpa_sub_dir <- file.path(base_output_dir, "VPA_Univariate", vpa_id)
  if (!dir.exists(vpa_sub_dir)) dir.create(vpa_sub_dir, recursive = TRUE)

  if (!is.null(vpa_conf$subset_to_groups)) {
    clin_df <- clin_df %>% filter(cohort_group %in% vpa_conf$subset_to_groups) %>% droplevels()
  }

  # EXCLUSION LOGGER & FILTERING
  all_cols <- setdiff(colnames(clin_df), "Subject_ID")
  vip_vars <- if(!is.null(vpa_conf$include_vars)) intersect(all_cols, vpa_conf$include_vars) else c()
  excl_log <- list()

  # 1. Exact Match Exclusions (No details needed)
  exact_excl <- setdiff(intersect(all_cols, vpa_conf$exclude_vars), vip_vars)
  if(length(exact_excl) > 0) {
    excl_log[["Exact Match (Blacklist)"]] <- exact_excl
    all_cols <- setdiff(all_cols, exact_excl)
  }

  # 2. Missingness, Variance, & Imbalance Exclusions (With Detailed Metrics)
  min_minority <- if(!is.null(vpa_conf$min_minority_pct)) vpa_conf$min_minority_pct else 0.10
  miss_excl <- c(); var_excl <- c(); imbalance_excl <- c(); surviving_cols <- c()

  for(v in all_cols) {
    if (v %in% vip_vars) {
      surviving_cols <- c(surviving_cols, v)
      next
    }

    vec <- clin_df[[v]]
    clean_vec <- vec[!is.na(vec)]
    miss_pct <- sum(is.na(vec)) / length(vec)

    # A. Missingness Filter
    if (miss_pct > vpa_conf$missingness_cutoff) {
      miss_excl <- c(miss_excl, sprintf("%s (%.1f%%)", v, miss_pct * 100))

      # B. Zero/Low Variance Filter
    } else if (length(unique(clean_vec)) < 2) {
      if (length(clean_vec) < 2) {
        var_excl <- c(var_excl, sprintf("%s (N < 2)", v))
      } else if (is.numeric(clean_vec)) {
        # Numeric variance is 0
        var_excl <- c(var_excl, sprintf("%s (Var: 0)", v))
      } else {
        # Categorical has only 1 level
        var_excl <- c(var_excl, sprintf("%s (<2 levels)", v))
      }

      # C. Extreme Imbalance Filter (Categorical only)
    } else if (is.factor(clean_vec) || is.character(clean_vec) || is.logical(clean_vec)) {
      props <- prop.table(table(clean_vec))
      min_prop <- min(props)
      min_level_name <- names(props)[which.min(props)]

      if (min_prop < min_minority) {
        imbalance_excl <- c(imbalance_excl, sprintf("%s (%s: %.1f%%)", v, min_level_name, min_prop * 100))
      } else {
        surviving_cols <- c(surviving_cols, v)
      }

      # Passed all filters
    } else {
      surviving_cols <- c(surviving_cols, v)
    }
  }

  if(length(miss_excl) > 0) excl_log[[sprintf("Missing > %s%%", vpa_conf$missingness_cutoff*100)]] <- miss_excl
  if(length(var_excl) > 0) excl_log[["Zero/Low Variance (<2 unique)"]] <- var_excl
  if(length(imbalance_excl) > 0) excl_log[[sprintf("Extreme Imbalance (< %s%% in minority class)", min_minority*100)]] <- imbalance_excl
  all_cols <- surviving_cols

  # 3. Partial Match Exclusions (No details needed)
  if(!is.null(vpa_conf$exclude_patterns)) {
    pat_excl <- c()
    for(pat in vpa_conf$exclude_patterns) {
      pat_excl <- c(pat_excl, grep(pat, all_cols, fixed = TRUE, value = TRUE))
    }
    pat_excl <- setdiff(unique(pat_excl), vip_vars)
    if(length(pat_excl) > 0) {
      excl_log[["Pattern Match (Substring)"]] <- pat_excl
      all_cols <- setdiff(all_cols, pat_excl)
    }
  }

  cols_to_scan <- all_cols
  if(length(cols_to_scan) == 0) return(NULL)

  excl_df <- data.frame(Reason = character(), Variables = character(), stringsAsFactors = FALSE)
  for(reason in names(excl_log)) {
    excl_df <- rbind(excl_df, data.frame(Reason = reason, Variables = paste(excl_log[[reason]], collapse = ", ")))
  }

  html_excl_log <- kableExtra::kbl(excl_df, format = "html", caption = "Variable Exclusion Log (with specific metrics)") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE) %>%
    kableExtra::column_spec(1, bold = TRUE, width = "25%")

  # UNIVARIATE SCAN (PARALLEL)
  results_list <- future_lapply(cols_to_scan, function(v) {
    vec <- clin_df[[v]]

    if (is.character(vec)) return(NULL)

    term <- if(is.factor(vec) || is.logical(vec)) {
      paste0("(1|`", v, "`)")
    } else {
      paste0("`", v, "`")
    }
    vpa_formula <- as.formula(paste("~", term))

    clin_tmp <- clin_df %>% select(Subject_ID, all_of(v)) %>% drop_na() %>% droplevels()
    valid_ids <- intersect(colnames(omics_mat), clin_tmp$Subject_ID)
    if(length(valid_ids) < 5) return(NULL)

    if (is.factor(clin_tmp[[v]]) && length(unique(clin_tmp[[v]])) < 2) return(NULL)
    if (is.numeric(clin_tmp[[v]]) && var(clin_tmp[[v]], na.rm = TRUE) == 0) return(NULL)

    res_var <- tryCatch({
      fitExtractVarPartModel(omics_mat[, valid_ids], vpa_formula, clin_tmp %>% filter(Subject_ID %in% valid_ids), BPPARAM = SerialParam())
    }, error = function(e) return(NULL))

    if(!is.null(res_var)) return(data.frame(Variable = v, Mean_Variance = mean(res_var[,1]), N_Samples = length(valid_ids)))
    return(NULL)
  }, future.seed = TRUE)

  full_results <- bind_rows(results_list) %>% arrange(desc(Mean_Variance))
  if (nrow(full_results) == 0) return(NULL)
  write_csv(full_results, file.path(vpa_sub_dir, paste0(vpa_id, "_Univariate_Results.csv")))

  # UNIVARIATE PLOT & SUMMARY TABLE
  plot_df <- head(full_results, vpa_conf$top_n_plot)

  p_scan <- ggplot(plot_df, aes(x = reorder(Variable, Mean_Variance), y = Mean_Variance)) +
    geom_bar(stat = "identity", fill = "#2c3e50", alpha = 0.8, color = "black") +
    geom_text(aes(label = sprintf("%.2f%%", Mean_Variance * 100)), hjust = -0.1, fontface = "bold") +
    coord_flip() + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +
    theme_project_base() + labs(title = vpa_conf$title, subtitle = "Univariate Leaderboard", y = "Mean Variance Explained", x = NULL)

  vars_to_summarize <- unique(plot_df$Variable)
  summary_rows <- lapply(vars_to_summarize, function(v) {
    if (!v %in% colnames(clin_df)) return(NULL)

    vec_full <- clin_df[[v]]
    miss_pct <- sum(is.na(vec_full)) / length(vec_full)
    vec <- vec_full[!is.na(vec_full)]

    var_val <- full_results$Mean_Variance[full_results$Variable == v]

    if (is.numeric(vec)) {
      details <- sprintf("Mean: %.2f | IQR: [%.2f - %.2f]", mean(vec), quantile(vec, 0.25), quantile(vec, 0.75))
      type <- "Continuous"
    } else {
      tbl <- prop.table(table(vec)) * 100
      lvl_strs <- paste0(names(tbl), ": ", sprintf("%.1f%%", tbl))
      if(length(lvl_strs) > 3) lvl_strs <- c(lvl_strs[1:3], "...(more)")
      details <- sprintf("Levels: %d | %s", length(table(vec)), paste(lvl_strs, collapse = ", "))
      type <- "Categorical"
    }

    return(data.frame(
      Variable = v,
      Type = type,
      Mean_Variance = scales::percent(var_val, accuracy = 0.01),
      Missingness = scales::percent(miss_pct, accuracy = 0.1),  # <-- NEW COLUMN HERE
      N_Valid = length(vec),
      Summary = details
    ))
  })

  html_summary_table <- bind_rows(summary_rows) %>%
    kableExtra::kbl(format = "html", caption = "Summary Statistics: Top Univariate Drivers") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE) %>%
    kableExtra::column_spec(3, bold = TRUE, color = "darkred")

  # RETURN PAYLOAD
  return_payload <- list(
    univariate_bar = p_scan,
    tables = list(exclusion = html_excl_log, summary = html_summary_table)
  )

  return(return_payload)
}

# 6. Pre-DEA Validation Wrapper & Formatter -----------------------------------
# Description:
#   Wraps the validate_dea_design mathematical engine. Parses the raw results
#   and formats them into a publication-ready, color-coded HTML diagnostic table.
# Args:
#   clin_df: Clinical metadata dataframe.
#   conf: A single DEA configuration list.
#   rules: The dea_validation_rules list containing thresholds.
#   cor_matrices: The loaded Spearman correlation matrices.
# Returns:
#   A list containing the raw engine results and the compiled kableExtra HTML table.
wrap_pre_dea_validation <- function(clin_df, conf, rules, cor_matrices = NULL) {

  res <- validate_dea_design(
    clin_df               = clin_df,
    group_col             = conf$group_col,
    target_groups         = conf$target_groups,
    ref_groups            = conf$ref_groups,
    covariates            = conf$covariates,
    min_samples_per_param = rules$min_samples_per_param,
    min_group_n           = rules$min_group_n,
    cor_matrices          = cor_matrices,
    max_cor_rho           = rules$max_cor_rho,
    max_cor_p             = rules$max_cor_p
  )

  status_color <- ifelse(res$Status == "Pass", "#27ae60", "#c0392b")

  spp_str <- case_when(
    is.na(res$SPP) ~ "-",
    res$SPP < rules$min_samples_per_param ~ sprintf("<span style='color:#c0392b; font-weight:bold;'>%.1f</span>", res$SPP),
    TRUE ~ as.character(round(res$SPP, 1))
  )

  params_str <- ifelse(is.na(res$Num_Params), "-", as.character(res$Num_Params))
  reasons_html <- paste(sprintf("&#8226; %s", res$Reasons), collapse = "<br>")

  matrix_rank_str <- case_when(
    all(res$Var_Details$Collinear == "-") ~ "<span style='color:gray;'>- (Skipped)</span>",
    any(grepl("Aliased", res$Var_Details$Collinear)) ~ "<span style='color:#c0392b; font-weight:bold;'>Rank Deficient (Aliased)</span>",
    TRUE ~ "<span style='color:black; font-weight:bold;'>Full Rank (Independent)</span>"
  )

  # FIX 4 & 1: Removed leading spaces so Markdown doesn't treat it as a code block
  meta_html <- sprintf(
    "<div style='text-align: left; padding: 15px; font-size: 14px; color: #333; background-color:#fcfcfc; border: 1px solid #ddd; margin-bottom: 15px;'>
<b>Cohort Survival:</b> Started with %d patients &rarr; <b>%d complete cases</b> survived.<br>
<b>Statistical Power:</b> <b>%s</b> samples per parameter (Total Parameters: %s).<br>
<b>Matrix Algebra:</b> %s<br><br>
<div style='margin-bottom: 5px;'><b>Verdict:</b> <span style='color:%s; font-size: 1.5em; font-weight: 900; letter-spacing: 1px;'>%s</span></div>
<b>Reason(s):</b><br>%s<br><br>
<i>Thresholds: Min Samples/Param = %d | Min Group N = %d | Max Corr = |%.2f|</i>
</div>\n\n",
    res$N_Start, res$N_Final, spp_str, params_str, matrix_rank_str,
    status_color, toupper(res$Status), reasons_html,
    rules$min_samples_per_param, rules$min_group_n, rules$max_cor_rho
  )

  plot_df <- res$Var_Details %>%
    mutate(
      Variance_Levels = case_when(
        grepl("^Var: 0\\.00$|^1 Levels$|^0 Levels$", Variance_Levels) ~ sprintf("<span style='color:#c0392b; font-weight:bold;'>%s</span>", Variance_Levels),
        TRUE ~ Variance_Levels
      ),
      High_Correlation = case_when(
        High_Correlation == "No" ~ "<span style='color:black;'>None</span>",
        High_Correlation == "-" ~ "<span style='color:gray;'>-</span>",
        TRUE ~ sprintf("<span style='color:#c0392b; font-weight:bold;'>%s</span>", High_Correlation)
      )
    ) %>%
    select(-Collinear)

  # FIX 4: Removed footnote from kable to render HTML separately
  tbl <- plot_df %>%
    kableExtra::kbl(
      format = "html",
      escape = FALSE,
      row.names = FALSE,
      col.names = c("Variable", "Role", "N (Valid)", "Missingness", "Distribution (Median [IQR] or %)", "High Correlation Warning"),
      caption = sprintf("<div style='color:black; font-size:1.4em; font-weight:bold; text-align:left; padding-bottom:5px;'>%s</div>", conf$title)
    ) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), full_width = TRUE, position = "left") %>%
    kableExtra::column_spec(1, bold = TRUE)

  # Passing meta_html as a separate object
  return(list(raw_data = res, html_table = tbl, meta_html = meta_html))
}

# 7. Wrapper for DEA Selector -------------------------------------------------
wrap_dea_selector <- function(clin_df, conf, vpa_df, cor_matrices) {

  res <- run_dea_selector(clin_df, conf, vpa_df, cor_matrices)
  if(res$status == "failed") return(list(status = "failed", msg = res$msg))

  spp_str <- if(res$meta$SPP < conf$min_samples_per_param) sprintf("<span style='color:#c0392b; font-weight:bold;'>%.1f (Overfitting Risk)</span>", res$meta$SPP) else sprintf("%.1f", res$meta$SPP)
  pct_cases <- (res$meta$N / res$meta$N_Start) * 100

  # FIX 1: Stripped leading spaces
  meta_html <- sprintf(
    "<div style='text-align: left; padding-top: 10px; font-size: 14px; color: #333;'>
<b>Final Complete Cases:</b> %d / %d (%.1f%%)<br>
<b>Total Parameters:</b> %d (Contrast + %d Covariates)<br>
<b>Statistical Power:</b> %s SPP (Threshold: %d)<br>
<b>Matrix Algebra:</b> %s<br>
<i>Algorithm Rules: Max Cor = %.2f, Max P = %.2f</i>
</div>",
    res$meta$N, res$meta$N_Start, pct_cases, res$meta$Params, length(res$selected_vars), spp_str, conf$min_samples_per_param, res$meta$Matrix, conf$max_cor_rho, conf$max_cor_p
  )

  formatted_tbl <- res$table %>%
    mutate(
      Recommendation = case_when(
        grepl("✅", Recommendation) ~ sprintf("<span style='color:#27ae60; font-weight:bold;'>%s</span>", Recommendation),
        TRUE ~ sprintf("<span style='color:#c0392b;'>%s</span>", Recommendation)
      ),
      Reason = ifelse(grepl("WARNING", Reason), sprintf("<span style='color:#e67e22; font-weight:bold;'>%s</span>", Reason), Reason)
    ) %>%
    kableExtra::kbl(format = "html", escape = FALSE, row.names = FALSE, caption = sprintf("<div style='color:black; font-size:1.4em; font-weight:bold;'>%s: Automated Selection</div>", conf$title)) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), full_width = TRUE) %>%
    kableExtra::footnote(general = meta_html, escape = FALSE, general_title = "")

  heatmap_vars <- c(conf$group_col, res$selected_vars)
  if (length(heatmap_vars) >= 2) {
    sub_cor <- cor_matrices$rho[heatmap_vars, heatmap_vars, drop = FALSE]
    sub_p <- cor_matrices$p[heatmap_vars, heatmap_vars, drop = FALSE]

    dist_mat <- as.dist(1 - sub_cor)
    hc <- hclust(dist_mat, method = "ward.D2")
    ordered_vars <- heatmap_vars[hc$order]

    plot_df <- as.data.frame(sub_cor) %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to="Var2", values_to="Correlation") %>%
      left_join(as.data.frame(sub_p) %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to="Var2", values_to="P_Value"), by=c("Var1", "Var2")) %>%
      mutate(
        idx_1 = match(Var1, ordered_vars), idx_2 = match(Var2, ordered_vars),
        Var2 = factor(Var2, levels = ordered_vars), Var1 = factor(Var1, levels = rev(ordered_vars)),
        Is_Sig = P_Value < conf$max_cor_p,
        Plot_Cor = if_else(Is_Sig, Correlation, NA_real_)
      ) %>% filter(idx_2 > idx_1)

    p_heatmap <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = Plot_Cor)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_x_discrete(position = "top") + scale_y_discrete(position = "right") +
      scale_fill_gradientn(colors = c(HM_Z_LOW, HM_Z_MID, HM_Z_HIGH), limits = c(-1, 1), na.value = "gray93", name = "Spearman Rho\n") +
      labs(title = sprintf("Correlation Heatmap: %s model", conf$title), x = NULL, y = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x.top = element_text(angle=45, hjust=0, vjust=0, margin=margin(b=5), face="bold"),
        axis.text.y.right = element_text(hjust=0, face="bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        plot.title.position = "plot"
      ) + coord_fixed()
  } else {
    p_heatmap <- NULL
  }

  return(list(
    status = "success",
    html_table = formatted_tbl,
    heatmap = p_heatmap,
    selected_vars = if(is.null(res$selected_vars)) character(0) else res$selected_vars
  ))
}

# 8. VPA Quality Control Pipeline Wrapper -------------------------------------
# Description:
#   Generates a distribution grid for the top VPA hits and a symmetrical
#   correlation heatmap sorted by VPA Rank to identify extreme skew and collinear blocks.
# Args:
#   vpa_df: The univariate VPA results dataframe.
#   clin_df: Clinical metadata dataframe.
#   cor_matrices: The loaded Spearman correlation matrices.
#   conf: The specific VPA configuration list.
#   vpa_id: String identifier for saving outputs.
#   base_output_dir: Path for saving results.
# Returns:
#   Nested list with 'plots' and 'html_warnings'.
wrap_vpa_qc_pipeline <- function(vpa_df, clin_df, cor_matrices, conf, vpa_id, base_output_dir) {
  results <- list(plots = list(), html_warnings = c())

  # Setup Directories
  qc_dir <- file.path(base_output_dir, "VPA_Univariate", "QC_Plots", vpa_id)
  if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

  clin_sub <- clin_df
  if (!is.null(conf$subset_to_groups)) {
    clin_sub <- clin_sub %>% filter(cohort_group %in% conf$subset_to_groups)
  }

  # A. Distribution Grid
  dist_vars <- head(vpa_df$Variable, conf$dist_plot_limit)
  dist_plots <- list()

  for (v in dist_vars) {
    if (!(v %in% colnames(clin_sub))) next
    vec <- clin_sub[[v]]

    if (is.numeric(vec)) {
      p <- ggplot(clin_sub, aes(x = .data[[v]])) +
        geom_histogram(fill = "grey50", color = "white", bins = 20, alpha = 0.8) +
        theme_project_base() + labs(title = v, x = NULL, y = "Count") +
        theme(plot.title = element_text(size = 10, face = "bold"))
    } else {
      p <- ggplot(clin_sub %>% filter(!is.na(.data[[v]])), aes(x = as.factor(.data[[v]]))) +
        geom_bar(fill = "grey70", color = "black", alpha = 0.8) +
        theme_project_base() + labs(title = v, x = NULL, y = "Count") +
        theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
    }
    dist_plots[[v]] <- p
  }

  if (length(dist_plots) > 0) {
    grid_plot <- wrap_plots(dist_plots, ncol = 4) +
      plot_annotation(title = sprintf("Distributions: Top %d VPA Hits (%s)", length(dist_plots), conf$title),
                      theme = theme(plot.title = element_text(size=16, face="bold")))
    results$plots[["distribution_grid"]] <- grid_plot

    ggsave(file.path(qc_dir, paste0(vpa_id, "_Distribution_Grid.png")), grid_plot, width = 12, height = 8, dpi = 300, bg="white")
  }

  # B. Correlation Heatmap
  # heat_vars_raw is naturally ordered by VPA Rank (highest first) + forced variables at the end
  heat_vars_raw <- unique(c(head(vpa_df$Variable, conf$corr_plot_limit), conf$heatmap_add_vars))
  available_vars <- colnames(cor_matrices$rho)

  missing_vars <- setdiff(heat_vars_raw, available_vars)

  # intersect() automatically preserves the order of the first argument (heat_vars_raw).
  # Therefore, valid_vars is already perfectly sorted by VPA Rank!
  valid_vars <- intersect(heat_vars_raw, available_vars)

  if (length(missing_vars) > 0) {
    results$html_warnings <- c(results$html_warnings, sprintf("> **WARNING:** Dropped from heatmap (Not strictly binarized/ordinal): *%s*", paste(missing_vars, collapse = ", ")))
  }

  if (length(valid_vars) >= 2) {
    sub_cor <- cor_matrices$rho[valid_vars, valid_vars, drop = FALSE]
    sub_p   <- cor_matrices$p[valid_vars, valid_vars, drop = FALSE]

    # BYPASS HIERARCHICAL CLUSTERING: Lock the visual order strictly to the VPA Rank
    ordered_vars <- valid_vars

    # Parameterized Threshold Fallback
    p_thresh <- if(exists("dea_validation_rules")) dea_validation_rules$max_cor_p else 0.05

    plot_df <- as.data.frame(sub_cor) %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to="Var2", values_to="Correlation") %>%
      left_join(as.data.frame(sub_p) %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to="Var2", values_to="P_Value"), by=c("Var1", "Var2")) %>%
      mutate(
        idx_1 = match(Var1, ordered_vars), idx_2 = match(Var2, ordered_vars),
        Var2 = factor(Var2, levels = ordered_vars), Var1 = factor(Var1, levels = rev(ordered_vars)),
        Is_Sig = P_Value < p_thresh,
        Plot_Cor = if_else(Is_Sig, Correlation, NA_real_)
      ) %>% filter(idx_2 > idx_1)

    p_heatmap <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = Plot_Cor)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_x_discrete(position = "top") + scale_y_discrete(position = "right") +
      scale_fill_gradientn(colors = c(HM_Z_LOW, HM_Z_MID, HM_Z_HIGH), limits = c(-1, 1), na.value = "gray93", name = "Spearman Rho\n") +
      # UI TWEAK: Update title to reflect the new sorting logic
      labs(title = sprintf("Correlation Network: Top %d Valid Hits (Rank Sorted)", length(valid_vars)), subtitle = sprintf("Grey squares indicate non-significance (p >= %s)", p_thresh), x = NULL, y = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x.top = element_text(angle=45, hjust=0, vjust=0, margin=margin(b=5), face="bold", size=rel(0.85)),
        axis.text.y.right = element_text(hjust=0, face="bold", size=rel(0.85)),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        plot.subtitle = element_text(color="grey40", hjust=0.5),
        plot.title.position = "plot"
      ) + coord_fixed()

    results$plots[["correlation_heatmap"]] <- p_heatmap

    # Dynamic sizing based on variable count
    dyn_dim <- max(6, length(valid_vars) * 0.4 + 2)
    ggsave(file.path(qc_dir, paste0(vpa_id, "_Correlation_Heatmap.png")), p_heatmap, width = dyn_dim, height = dyn_dim, dpi = 300, bg="white")
  }

  return(results)
}

# 9. VPA Scatter Diagnostics Engine -------------------------------------------
# Description:
#   Extracts the Top 6 proteins driven by each covariate and plots their NPX
#   distributions (scatter + regression for numeric, boxplots for factors) to
#   visually check for leverage points and outliers.
wrap_vpa_scatter_diagnostics <- function(omics_mat, clin_df, vpa_df, covariates, model_id, title_prefix, base_output_dir) {
  results <- list(plots = list(), paths = list())

  diag_dir <- file.path(base_output_dir, "VPA_Diagnostics", model_id)
  if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

  for (cov in covariates) {
    if (!(cov %in% colnames(clin_df)) || !(cov %in% colnames(vpa_df))) next

    # Get Top 6 Proteins for THIS specific covariate
    top_feats <- vpa_df %>% arrange(desc(.data[[cov]])) %>% head(6) %>% pull(Feature)
    if(length(top_feats) == 0) next

    # Extract & Filter Clinical Data
    plot_data <- clin_df %>% select(Subject_ID, cov_val = all_of(cov)) %>% drop_na()
    valid_ids <- intersect(plot_data$Subject_ID, colnames(omics_mat))
    plot_data <- plot_data %>% filter(Subject_ID %in% valid_ids)

    # Extract & Melt Omics Data
    omics_sub <- omics_mat[top_feats, valid_ids, drop = FALSE]
    omics_long <- as.data.frame(t(omics_sub)) %>% rownames_to_column("Subject_ID") %>%
      pivot_longer(-Subject_ID, names_to = "Protein", values_to = "NPX")

    plot_df <- plot_data %>% left_join(omics_long, by = "Subject_ID") %>%
      mutate(Protein = factor(Protein, levels = top_feats)) # Preserve VPA rank order

    safe_cov_name <- gsub("[^A-Za-z0-9_.-]", "_", cov)

    # Generate Plot based on Data Type
    if (is.numeric(plot_df$cov_val)) {
      p <- ggplot(plot_df, aes(x = cov_val, y = NPX, color = Protein)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
        # Add R-Squared text to each facet
        ggpubr::stat_cor(aes(label = after_stat(rr.label)), color = "black", size = 3.5, label.y.npc = "top") +
        facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
        theme_project_base() +
        theme(legend.position = "none") +
        labs(title = sprintf("%s", title_prefix),
             subtitle = sprintf("Top 6 Variance Drivers for '%s' (Free Y-Axis)", cov),
             x = cov, y = "NPX Value")
    } else {
      p <- ggplot(plot_df, aes(x = as.factor(cov_val), y = NPX, fill = Protein)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        geom_jitter(aes(color = Protein), width = 0.2, alpha = 0.7, size = 1.5) +
        facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
        theme_project_base() +
        theme(legend.position = "none") +
        labs(title = sprintf("%s", title_prefix),
             subtitle = sprintf("Top 6 Variance Drivers for '%s' (Free Y-Axis)", cov),
             x = cov, y = "NPX Value")
    }

    plot_path <- file.path(diag_dir, paste0(model_id, "_", safe_cov_name, ".png"))
    ggsave(plot_path, p, width = 11, height = 7, dpi = 300, bg = "white")

    results$plots[[cov]] <- p
    results$paths[[cov]] <- plot_path
  }

  return(results)
}




# 10. Matrix Factory Engine ----------------------------------------------------
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

# 12. Consensus Export & Alignment Engine --------------------------------------
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


# 13. Global Omics Variance Filter Engine --------------------------------------
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

# 14. Clinical Phenotyping Engine (ComplexHeatmap version) ---------------------
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

# 15. K-Optimization (Pre-Flight) Pipeline Wrapper ---------------------------
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

# 16. UMAP & Multi-Algorithm Molecular Clustering Engine -----------------------
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


# 17. MISSING DATA IMPUTATION: RANDOM FOREST (missForest) ----------------------
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


# 18. GOWER DISTANCE FACTORY (For Mixed Data) ----------------------------------
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

# 19. FAMD PIPELINE (Factor Analysis of Mixed Data) ----------------------------
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

# PART 2: INTERNAL HELPER ENGINES (The "Workers") ------------------------------

# 1. Limma & Volcano Engine ----------------------------------------------------
# Description:
#   The mathematical core for differential expression. Creates a design matrix
#   (including any specified covariates), fits linear models using eBayes, and
#   constructs a Volcano plot annotated with the top changed features.
# Args:
#   omics_mat: Feature matrix.
#   clin_df: Metadata (already filtered by wrapper).
#   dea_conf: DEA configuration lists.
#   project_colors: Global color function for plotting.
# Returns:
#   List containing the topTable dataframe and the ggplot Volcano object.
run_limma_engine <- function(omics_mat, clin_df, dea_conf, project_colors) {

  required_cols <- c(dea_conf$group_col, dea_conf$covariates)
  clin_sub <- clin_df %>% drop_na(all_of(required_cols))

  clin_sub <- clin_sub %>% mutate(limma_group = case_when(
    .data[[dea_conf$group_col]] %in% dea_conf$target_groups ~ "Target",
    .data[[dea_conf$group_col]] %in% dea_conf$ref_groups ~ "Reference",
    TRUE ~ NA_character_
  )) %>%
    filter(!is.na(limma_group)) %>%
    mutate(limma_group = factor(limma_group, levels = c("Reference", "Target")))

  if(length(unique(clin_sub$limma_group)) < 2) return(NULL)

  mat_sub <- omics_mat[, clin_sub$Subject_ID, drop = FALSE]

  if (length(dea_conf$covariates) > 0) {
    # Wraps all covariates in backticks to protect special characters like -, %, ()
    protected_covariates <- paste(sprintf("`%s`", dea_conf$covariates), collapse = " + ")
    formula_str <- paste("~ 0 + limma_group +", protected_covariates)
  } else {
    formula_str <- "~ 0 + limma_group"
  }

  design <- model.matrix(as.formula(formula_str), data = clin_sub)
  colnames(design)[1:2] <- c("Reference", "Target")

  # --- Sanitize all column names so limma's strict parser doesn't crash ---
  colnames(design) <- make.names(colnames(design))

  fit <- lmFit(mat_sub, design)
  contrast_mat <- makeContrasts(Target - Reference, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_mat))

  dea_res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P") %>%
    rownames_to_column("Feature") %>%
    mutate(
      Significance = case_when(
        adj.P.Val < dea_conf$dea_p_cutoff & logFC > dea_conf$dea_fc_cutoff  ~ "Up-regulated",
        adj.P.Val < dea_conf$dea_p_cutoff & logFC < -dea_conf$dea_fc_cutoff ~ "Down-regulated",
        TRUE ~ "Not Significant"
      ),
      label_score = abs(logFC) * -log10(adj.P.Val)
    )

  top_labels <- dea_res %>%
    filter(Significance != "Not Significant") %>%
    slice_max(label_score, n = 20, with_ties = FALSE)

  target_color <- project_colors(dea_conf$target_groups[1])
  ref_color <- project_colors(dea_conf$ref_groups[1])

  p_volc <- ggplot(dea_res, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down-regulated" = HM_Z_LOW, "Not Significant" = "grey85", "Up-regulated" = HM_Z_HIGH)) +
    geom_text_repel(data = top_labels, aes(label = Feature), size = 3, color = "black", box.padding = 0.5, max.overlaps = Inf) +
    geom_vline(xintercept = c(-dea_conf$dea_fc_cutoff, dea_conf$dea_fc_cutoff), linetype = "dotted") +
    geom_hline(yintercept = -log10(dea_conf$dea_p_cutoff), linetype = "dashed", color = "darkred") +
    theme_project_base() +
    labs(title = sprintf("Volcano: %s", dea_conf$title), x = "Log2 Fold Change", y = "-log10(FDR)")

  return(list(data = dea_res, volcano = p_volc))
}

# 2. GSEA Engine ---------------------------------------------------------------
# Description:
#   Ranks the Limma features by t-statistic and performs Gene Set Enrichment
#   Analysis. Automatically cleans pathway names and generates DotPlots and RidgePlots.
# Args:
#   dea_res: topTable dataframe from run_limma_engine.
#   go_db: Loaded database (e.g., MSigDB or generic TERM2GENE table).
#   omics_config: Settings (p-value, min sizes).
#   title: Contrast name for plot headers.
# Returns:
#   Dataframe of pathways, a dotplot ggplot, and a ridgeplot ggplot.
run_gsea_engine <- function(dea_res, go_db, omics_config, title, db_name) {

  if (is.null(go_db)) return(NULL)

  ranked_vec <- setNames(dea_res$t, dea_res$Feature) %>% sort(decreasing = TRUE)
  set.seed(42)

  # Calculate GSEA for ALL pathways
  gsea_res <- clusterProfiler::GSEA(geneList = ranked_vec, TERM2GENE = go_db, pvalueCutoff = 1, minGSSize = omics_config$gsea_min_size, verbose = FALSE)
  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(NULL)

  gsea_res@result$Description <- gsub("GOBP_", "", gsea_res@result$Description) %>% gsub("_", " ", .)
  res_df <- as.data.frame(gsea_res) %>% mutate(Status = ifelse(NES > 0, "Up-regulated", "Down-regulated"))

  warning_tag <- if(sum(res_df$p.adjust < omics_config$gsea_p_cutoff) == 0) "\n(Exploratory: No paths passed FDR)" else ""

  # Wrap long geneset names for the Dotplot
  plot_df <- res_df %>%
    group_by(Status) %>%
    slice_max(abs(NES), n = 10) %>%
    ungroup() %>%
    mutate(Description = stringr::str_wrap(Description, width = 45))

  main_title <- str_wrap(sprintf("GSEA (%s): %s", db_name, title), width = 60)
  full_title <- paste0(main_title, warning_tag)

  # 1. Exploratory Dotplot
  p_dot <- ggplot(plot_df, aes(x = NES, y = reorder(Description, NES), color = p.adjust, size = setSize)) +
    geom_point() +
    scale_color_gradient(low = "firebrick3", high = "navy") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme_project_base() +
    theme(plot.title.position = "plot", plot.title = element_text(hjust = 0)) +
    labs(title = full_title, y = NULL)

  # STRICT BIOLOGICAL PRE-FILTERING (FDR < cutoff) for Networks
  gsea_sig <- gsea_res
  gsea_sig@result <- gsea_sig@result %>% filter(p.adjust < omics_config$network_p_cutoff)

  # NOTE: simplify() has been permanently bypassed here due to generic GSEA() incompatibility.

  gsea_sim <- if(nrow(gsea_sig@result) >= 2) enrichplot::pairwise_termsim(gsea_sig) else NULL

  # 2. Ridgeplot
  ridge_title <- str_wrap(sprintf("Pathways (%s): %s | FDR cutoff: %.2f", db_name, title, omics_config$network_p_cutoff), width = 60)
  p_ridge <- if(nrow(gsea_sig@result) > 0) {
    enrichplot::ridgeplot(gsea_sig, showCategory = 15) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
      theme_project_base() +
      labs(title = ridge_title)
  } else NULL

  # 3. Enrichment Map (MATH SHIELD ACTIVATED)
  p_emap <- NULL
  if(!is.null(gsea_sim)) {
    sim_mat <- gsea_sim@termsim
    if (!is.null(sim_mat) && nrow(sim_mat) > 1) {
      valid_edges <- sum(sim_mat[upper.tri(sim_mat)] >= 0.3, na.rm = TRUE)

      if (valid_edges > 0) {
        emap_title <- str_wrap(sprintf("Enrichment Map (%s): %s | FDR cutoff: %.2f", db_name, title, omics_config$network_p_cutoff), width = 60)
        p_emap <- enrichplot::emapplot(
          gsea_sim, color = "NES", showCategory = 40, node_label_size = 2.5,
          size_category = 0.8, size_edge = 0.1, color_edge = "grey60", min_edge = 0.3
        ) +
          theme_project_base() +
          labs(title = emap_title, subtitle = "Nodes = Pathways | Edges = Shared Genes")
      }
    }
  }

  return(list(data = res_df, dotplot = p_dot, ridgeplot = p_ridge, emap = p_emap))
}

# 3. Heatmap Engine ------------------------------------------------------------
# Description:
#   Selects the top N most significant features from a DEA result, converts them
#   to Z-scores across the requested subjects, and renders an annotated heatmap
#   using the ComplexHeatmap library.
# Args:
#   mat: Raw omics matrix.
#   clin: Metadata for annotations.
#   dea_res: Dataframe of Limma results used to select features.
#   title: Plot title.
#   n_top: Number of rows to include in the heatmap.
#   split_by_group: Boolean dictating if the columns should be visually clustered by cohort.
# Returns:
#   A ComplexHeatmap object (must be drawn with draw()).
run_heatmap_engine <- function(mat, clin, dea_res, title, n_top = 50, split_by_group = FALSE, split_col = "cohort_group", show_top_n_title = TRUE, custom_color_map = NULL, annotation_vars = NULL) {

  valid_ids <- intersect(clin$Subject_ID, colnames(mat))
  clin_sub <- clin %>% filter(Subject_ID %in% valid_ids)

  top_feats <- dea_res %>% mutate(pi = abs(logFC) * -log10(adj.P.Val)) %>% slice_max(pi, n = n_top) %>% pull(Feature)

  if (length(top_feats) < 2) return(NULL) # Safety catch

  plot_mat <- t(scale(t(mat[top_feats, clin_sub$Subject_ID])))


  # 1. Determine which variables to extract for the top annotation
  if (is.null(annotation_vars)) {
    # Legacy Support: Default to the split column, Sex, and Age
    req_vars <- unique(c(split_col, "Sex", "Age"))
  } else {
    # Dynamic Mode: Ensure the split column is always included to visualize the contrast
    req_vars <- unique(c(split_col, annotation_vars))
  }

  valid_vars <- intersect(req_vars, colnames(clin_sub))

  # 2. Build the Annotation DataFrame
  anno_df <- clin_sub %>%
    select(all_of(valid_vars)) %>%
    mutate(across(where(is.character), as.factor))

  # 3. Dynamically Assign Colors based on your global constants
  col_list <- list()

  # Handle the primary split column (e.g., Dynamic_Cluster)
  if (!is.null(split_col) && split_col %in% valid_vars) {
    col_list[[split_col]] <- get_project_colors(as.character(unique(clin_sub[[split_col]])), custom_map = custom_color_map)
  }

  # EXPLICIT TRAP: Force cohort_group to ALWAYS use the official project palette
  if ("cohort_group" %in% valid_vars) {
    col_list[["cohort_group"]] <- get_project_colors(as.character(unique(clin_sub[["cohort_group"]])))
  }

  # Explicit Trap: Sex (Now correctly utilizing get_project_colors)
  if ("Sex" %in% valid_vars) {
    col_list[["Sex"]] <- get_project_colors(as.character(unique(clin_sub[["Sex"]])))
  }

  # Explicit Trap: Age
  if ("Age" %in% valid_vars) {
    col_list[["Age"]] <- circlize::colorRamp2(c(min(clin_sub$Age, na.rm=T), max(clin_sub$Age, na.rm=T)), c(COLOR_AGE_LOW, COLOR_AGE_HIGH))
  }



  # 4. Construct the HeatmapAnnotation object dynamically
  if (ncol(anno_df) > 0) {
    ha <- do.call(ComplexHeatmap::HeatmapAnnotation,
                  c(as.list(anno_df),
                    list(col = if(length(col_list) > 0) col_list else NULL,
                         na_col = "grey80",
                         annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"))))
  } else {
    ha <- NULL
  }

  final_title <- if(show_top_n_title) sprintf("%s: Top %d Proteins", title, n_top) else title

  # 5. Draw the Heatmap
  ComplexHeatmap::Heatmap(plot_mat, name = "Z-Score", column_title = final_title,
                          top_annotation = ha, show_column_names = FALSE, column_split = if(split_by_group) clin_sub[[split_col]] else NULL,
                          col = circlize::colorRamp2(c(-2, 0, 2), c(HM_Z_LOW, HM_Z_MID, HM_Z_HIGH)))
}

# 4. LASSO Engine --------------------------------------------------------------
# Description:
#   Runs a generalized linear model with L1 penalty (LASSO) using cross-validation.
#   Extracts the coefficients at the optimal lambda ("lambda.min") to identify
#   which features best discriminate between the two target classes.
# Args:
#   mat: Omics matrix.
#   clin: Clinical dataframe.
#   target_col: Name of the column containing class labels.
#   target_groups: The two labels to classify (e.g., c("Active", "Control")).
#   title: Plot title prefix.
# Returns:
#   List containing the glmnet cross-validation object, a bar chart of non-zero
#   coefficients, and the coefficient dataframe.
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

# 5. Cross-Contrast Engine ----------------------------------------------------
# Description:
#   Joins two different Limma topTable dataframes on 'Feature'. Determines if
#   features are significant in one, the other, or both. Maps them to a quadrant
#   scatter plot to reveal biological divergence or shared mechanisms.
# Args:
#   dea_x, dea_y: Limma result dataframes.
#   conf: Configuration lists with 'top_n' parameters for labeling.
#   x_label, y_label: Axis and legend names.
#   color_x, color_y: Highlight colors for specific signatures.
#   highlight_color: Color for shared signatures.
#   title: Plot title.
# Returns:
#   A fully formatted 4-quadrant ggplot scatter.
run_cross_contrast_engine <- function(dea_x, dea_y, conf, x_label, y_label, color_x, color_y, highlight_color, title) {
  x_sig_name <- sprintf("%s Specific", x_label)
  y_sig_name <- sprintf("%s Specific", y_label)

  temporal_shifts <- full_join(
    dea_x %>% select(Feature, logFC_X = logFC, FDR_X = adj.P.Val, Sig_Status_X = Significance),
    dea_y %>% select(Feature, logFC_Y = logFC, FDR_Y = adj.P.Val, Sig_Status_Y = Significance),
    by = "Feature"
  ) %>%
    mutate(
      Sig_X = !is.na(Sig_Status_X) & Sig_Status_X != "Not Significant",
      Sig_Y = !is.na(Sig_Status_Y) & Sig_Status_Y != "Not Significant",
      Signature = case_when(
        Sig_X & Sig_Y  ~ "Shared Drivers",
        Sig_X & !Sig_Y ~ x_sig_name,
        !Sig_X & Sig_Y ~ y_sig_name,
        TRUE ~ "Not Significant"
      ),
      label_score = case_when(
        Signature == "Shared Drivers" ~ (abs(logFC_X) * -log10(FDR_X)) + (abs(logFC_Y) * -log10(FDR_Y)),
        Signature == x_sig_name ~ abs(logFC_X) * -log10(FDR_X),
        Signature == y_sig_name ~ abs(logFC_Y) * -log10(FDR_Y),
        TRUE ~ 0
      ),
      Quadrant = case_when(
        logFC_X > 0 & logFC_Y > 0 ~ "Up-Up",
        logFC_X < 0 & logFC_Y < 0 ~ "Down-Down",
        logFC_X > 0 & logFC_Y < 0 ~ "Up-Down",
        logFC_X < 0 & logFC_Y > 0 ~ "Down-Up",
        TRUE ~ "None"
      )
    )

  sig_points <- temporal_shifts %>% filter(Signature != "Not Significant")
  if(nrow(sig_points) == 0) return(NULL)

  top_n_spec <- if(!is.null(conf$top_n_specific)) conf$top_n_specific else 5
  top_n_shar <- if(!is.null(conf$top_n_shared)) conf$top_n_shared else 5

  labels_shared <- sig_points %>% filter(Signature == "Shared Drivers") %>% group_by(Quadrant) %>% slice_max(order_by = label_score, n = top_n_shar, with_ties = FALSE) %>% ungroup()
  labels_specific <- sig_points %>% filter(Signature %in% c(x_sig_name, y_sig_name)) %>% group_by(Signature) %>% slice_max(order_by = label_score, n = top_n_spec, with_ties = FALSE) %>% ungroup()
  top_labels <- bind_rows(labels_shared, labels_specific)

  sig_points <- sig_points %>% mutate(repel_label = ifelse(Feature %in% top_labels$Feature, Feature, ""))

  color_mapping <- setNames(c(highlight_color, color_x, color_y), c("Shared Drivers", x_sig_name, y_sig_name))

  ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = 0, color = "grey40") + geom_vline(xintercept = 0, color = "grey40") +
    geom_point(data = sig_points, aes(x = logFC_X, y = logFC_Y, color = Signature), alpha = 0.8, size = 2.5) +
    geom_text_repel(data = sig_points, aes(x = logFC_X, y = logFC_Y, label = repel_label), color = "black", size = 3.5, box.padding = 0.8, max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = color_mapping) +
    theme_project_base() +
    labs(title = title, x = paste(x_label, "(Log2 Fold Change)"), y = paste(y_label, "(Log2 Fold Change)")) +
    coord_fixed(ratio = 1)
}

# 6. Shared VPA Worker Engine -------------------------------------------------
# Description:
#   The math core for Variance Partitioning. Uses linear mixed models via the
#   variancePartition package. Calculates the mean variance and 95% CI explained
#   by factors defined in the formula.
# Args:
#   mat: Matrix of subset omics data.
#   clin: Dataframe of aligned metadata.
#   vpa_formula: Fully constructed R formula string (e.g., ~ Age + (1|Sex)).
#   title: Plot title.
#   top_n: Number of variables to display in the plot limit.
# Returns:
#   Raw variance dataframe, a barplot (with CI bounds), and a violin distribution plot.
run_vpa_engine <- function(mat, clin, vpa_formula, title, top_n = 10) {

  # Fit the Model (PARALLELIZED)
  varPart <- tryCatch({
    fitExtractVarPartModel(mat, vpa_formula, clin, BPPARAM = SnowParam(workers = 4))
  }, error = function(e) return(NULL))

  if(is.null(varPart)) return(NULL)
  varPart_df <- as.data.frame(varPart) %>% rownames_to_column("Feature")

  # Calculate Statistics
  summary_df <- data.frame(
    Covariate = colnames(varPart),
    Mean_Var = colMeans(varPart, na.rm = TRUE),
    SEM = apply(varPart, 2, sd, na.rm = TRUE) / sqrt(nrow(varPart))
  ) %>%
    mutate(CI_upper = pmin(1, Mean_Var + (1.96 * SEM))) %>%
    arrange(desc(Mean_Var))

  plot_summary <- head(summary_df, top_n)

  # Barplot
  p_bar <- ggplot(plot_summary, aes(x = reorder(Covariate, Mean_Var), y = Mean_Var, fill = Covariate)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, Mean_Var - (1.96 * SEM)), ymax = CI_upper), width = 0.2) +
    geom_text(aes(y = CI_upper, label = sprintf("%.1f%%", Mean_Var * 100)), hjust = -0.2, fontface = "bold") +
    coord_flip() + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +
    theme_project_base() + theme(legend.position = "none") +
    labs(title = title, subtitle = "Average Variance Explained", y = "% Variance Explained (± 95% CI)", x = NULL)

  # Violin Plot
  long_var <- varPart_df %>%
    select(Feature, all_of(plot_summary$Covariate)) %>%
    pivot_longer(-Feature, names_to = "Covariate", values_to = "Val") %>%
    mutate(Covariate = factor(Covariate, levels = rev(plot_summary$Covariate)))

  p_violin <- ggplot(long_var, aes(x = Covariate, y = Val, fill = Covariate)) +
    geom_violin(alpha = 0.6, color = "grey40", scale = "width") +
    geom_boxplot(width = 0.08, fill = "white", color = "black", outlier.shape = 19, outlier.size = 1, outlier.alpha = 0.4) +
    coord_flip() + scale_y_continuous(labels = scales::percent) +
    theme_project_base() + theme(legend.position = "none") +
    labs(title = title, x = NULL, y = "% Variance Explained")

  return(list(data = varPart_df, plots = list(bar = p_bar, violin = p_violin)))
}

# 7. Pre-DEA Validation Engine (Variable-Level Profiler) ----------------------
validate_dea_design <- function(clin_df, group_col, target_groups, ref_groups, covariates,
                                min_samples_per_param = 10, min_group_n = 2,
                                cor_matrices = NULL, max_cor_rho = 0.7, max_cor_p = 0.05) {

  all_vars <- c(group_col, covariates)

  var_details <- data.frame(
    Variable = all_vars,
    Role = c("Contrast Group", rep("Covariate", length(covariates))),
    N_Valid = NA_integer_,
    Missingness = NA_character_,
    Variance_Levels = NA_character_,
    Collinear = "-",
    High_Correlation = "-",
    stringsAsFactors = FALSE
  )

  # NEW SAFETY CATCH: Check if variables actually exist in the dataframe
  missing_vars <- setdiff(all_vars, colnames(clin_df))
  if (length(missing_vars) > 0) {
    return(list(Status = "Fail",
                Reasons = c(sprintf("FATAL: Variables not found in clinical data: %s", paste(missing_vars, collapse = ", "))),
                Var_Details = var_details, N_Start = 0, N_Final = 0, SPP = NA, Num_Params = NA))
  }

  sub_df <- clin_df %>% filter(!!sym(group_col) %in% c(target_groups, ref_groups))
  n_start <- nrow(sub_df)

  if (n_start == 0) return(list(Status = "Fail", Reasons = c("Target/Ref groups not found in data."), Var_Details = var_details, N_Start = 0, N_Final = 0, SPP = NA, Num_Params = NA))

  # 2. Variable-Level N_Valid & Missingness
  for(i in 1:nrow(var_details)) {
    v <- var_details$Variable[i]
    valid_count <- sum(!is.na(sub_df[[v]]))
    var_details$N_Valid[i] <- valid_count

    pct <- (n_start - valid_count) / n_start
    var_details$Missingness[i] <- sprintf("%.1f%%", pct * 100)
  }

  # 3. Complete Case Filter
  cc_df <- sub_df %>% select(Subject_ID, all_of(all_vars)) %>% drop_na() %>% droplevels()
  n_final <- nrow(cc_df)

  if (n_final == 0) return(list(Status = "Fail", Reasons = c("100% NA-Wipeout. No patients have complete data."), Var_Details = var_details, N_Start = n_start, N_Final = 0, SPP = NA, Num_Params = NA))

  reasons <- c()
  has_zero_variance <- FALSE

  # 4. Variable-Level Variance & Group Limits (UPDATED TO DESCRIPTIVE STATS)
  for(i in 1:nrow(var_details)) {
    v <- var_details$Variable[i]
    vec <- cc_df[[v]]

    if (is.numeric(vec)) {
      var_val <- var(vec, na.rm = TRUE)
      if (var_val == 0) {
        var_details$Variance_Levels[i] <- "Var: 0.00"
        reasons <- c(reasons, sprintf("'%s' has 0 variance.", v))
        has_zero_variance <- TRUE
      } else {
        # Valid continuous variable: show Median and IQR
        med_val <- median(vec, na.rm = TRUE)
        q25 <- quantile(vec, 0.25, na.rm = TRUE)
        q75 <- quantile(vec, 0.75, na.rm = TRUE)
        var_details$Variance_Levels[i] <- sprintf("Med: %.1f [%.1f - %.1f]", med_val, q25, q75)
      }
    } else {
      lvl_count <- length(unique(vec))
      if (lvl_count < 2) {
        var_details$Variance_Levels[i] <- sprintf("%d Levels", lvl_count)
        reasons <- c(reasons, sprintf("'%s' collapsed to <2 levels.", v))
        has_zero_variance <- TRUE
      } else {
        # Valid categorical variable: show Breakdown percentages
        tbl <- prop.table(table(vec)) * 100
        lvl_strs <- paste0(names(tbl), ": ", sprintf("%.1f%%", tbl))
        var_details$Variance_Levels[i] <- paste(lvl_strs, collapse = "<br>")
      }
    }
  }

  t_count <- sum(cc_df[[group_col]] %in% target_groups)
  r_count <- sum(cc_df[[group_col]] %in% ref_groups)
  if (t_count < min_group_n || r_count < min_group_n) {
    reasons <- c(reasons, sprintf("Primary groups too small (Target: %d, Ref: %d. Min: %d).", t_count, r_count, min_group_n))
    has_zero_variance <- TRUE
  }

  num_params <- NA
  spp <- NA

  # 5. Collinearity, Correlation, & Overfitting
  if (!has_zero_variance) {
    var_details$Collinear <- "No"
    var_details$High_Correlation <- "No"

    # A. Matrix Rank Check (Fatal Collinearity)
    cc_df[[group_col]] <- factor(cc_df[[group_col]], levels = c(ref_groups[1], target_groups[1]))
    formula_str <- paste("~", paste(sprintf("`%s`", all_vars), collapse = " + "))

    design <- tryCatch(model.matrix(as.formula(formula_str), data = cc_df), error = function(e) NULL)

    if (is.null(design)) {
      reasons <- c(reasons, "Could not construct linear design matrix.")
    } else {
      num_params <- ncol(design)
      q <- qr(design)

      if (q$rank < num_params) {
        aliased_cols <- colnames(design)[q$pivot[(q$rank + 1):num_params]]
        aliased_vars <- c()

        for(i in 1:nrow(var_details)) {
          v <- var_details$Variable[i]
          if (any(grepl(v, aliased_cols, fixed = TRUE))) {
            var_details$Collinear[i] <- "Yes (Aliased)"
            aliased_vars <- c(aliased_vars, v)
          }
        }
        reasons <- c(reasons, sprintf("Algebraic Collinearity detected. Aliased variable(s): %s.", paste(aliased_vars, collapse = ", ")))
      }

      # B. Spearman Variance Inflation Check (Biological Correlation) - UPDATED TO INCLUDE EXACT VALUES
      if (!is.null(cor_matrices) && length(covariates) > 1) {
        cor_vars <- intersect(covariates, rownames(cor_matrices$rho))
        if (length(cor_vars) > 1) {
          combos <- combn(cor_vars, 2, simplify = FALSE)
          for (combo in combos) {
            v1 <- combo[1]; v2 <- combo[2]
            rho_val <- cor_matrices$rho[v1, v2]
            p_val <- cor_matrices$p[v1, v2]

            if (!is.na(rho_val) && !is.na(p_val)) {
              if (abs(rho_val) >= max_cor_rho && p_val <= max_cor_p) {
                reasons <- c(reasons, sprintf("Variance Inflation Risk: '%s' & '%s' are highly correlated (rho = %.2f).", v1, v2, rho_val))

                # Append exactly the details to the output column
                exact_warning <- sprintf("Yes (with %s, &rho; = %.2f, p = %.3f)", v2, rho_val, p_val)
                exact_warning2 <- sprintf("Yes (with %s, &rho; = %.2f, p = %.3f)", v1, rho_val, p_val)

                var_details$High_Correlation[var_details$Variable == v1] <- exact_warning
                var_details$High_Correlation[var_details$Variable == v2] <- exact_warning2
              }
            }
          }
        }
      }

      # C. Overfitting Check
      spp <- n_final / num_params
      if (spp < min_samples_per_param) {
        reasons <- c(reasons, sprintf("Overfitting risk. Only %.1f samples per parameter (Requires >= %d).", spp, min_samples_per_param))
      }
    }
  } else {
    reasons <- c(reasons, "Skipped collinearity & correlation checks due to missing groups or zero variance.")
  }

  status <- if(length(reasons) == 0) "Pass" else "Fail"
  if (length(reasons) == 0) reasons <- c("All metrics within acceptable mathematical thresholds.")

  return(list(Status = status, Reasons = reasons, Var_Details = var_details, N_Start = n_start, N_Final = n_final, SPP = spp, Num_Params = num_params))
}

# 8. Automated DEA Variable Selector (Greedy Algorithm) -----------------------
run_dea_selector <- function(clin_df, conf, vpa_df, cor_matrices) {

  group_col <- conf$group_col
  sub_df <- clin_df %>% filter(!!sym(group_col) %in% c(conf$target_groups, conf$ref_groups))
  n_start <- nrow(sub_df)

  if (n_start == 0) return(list(status = "failed", msg = "Target/Ref groups not found."))

  # Ensure target column is numeric for correlation
  sub_df[[group_col]] <- as.numeric(factor(sub_df[[group_col]], levels = c(conf$ref_groups[1], conf$target_groups[1]))) - 1

  # Prep Candidate List: must_include first, then remaining sorted by VPA
  candidates <- unique(c(conf$must_include, vpa_df$Variable))
  candidates <- setdiff(candidates, c(group_col, conf$must_exclude))

  results_tbl <- data.frame(
    Variable = character(), VPA_Rank = character(), Distribution = character(),
    Recommendation = character(), Reason = character(), VIF = character(),
    stringsAsFactors = FALSE
  )

  selected_covariates <- c()
  current_cc_n <- n_start

  for (i in seq_along(candidates)) {
    v <- candidates[i]
    is_forced <- v %in% conf$must_include
    vpa_rank <- ifelse(v %in% vpa_df$Variable, sprintf("#%d", which(vpa_df$Variable == v)), "-")

    vec <- sub_df[[v]]

    # UI Tweak: Hide massive factor distributions
    dist_str <- if(is.numeric(vec)) {
      sprintf("%.1f [%.1f - %.1f]", median(vec, na.rm=T), quantile(vec, 0.25, na.rm=T), quantile(vec, 0.75, na.rm=T))
    } else {
      tbl <- prop.table(table(vec)) * 100
      if(length(tbl) > 4) {
        "<span style='color:grey;'>Too many levels to display</span>"
      } else {
        paste0(names(tbl), ": ", sprintf("%.1f%%", tbl), collapse = "<br>")
      }
    }

    # 1. Complete Case Check
    temp_vars <- c(group_col, selected_covariates, v)
    temp_n <- sub_df %>% select(Subject_ID, all_of(temp_vars)) %>% drop_na() %>% nrow()

    # 2. Power Check (Rule of 10) - UI Tweak: General Early Stop Message
    spp <- temp_n / length(temp_vars)
    if (!is_forced && spp < conf$min_samples_per_param) {
      results_tbl[nrow(results_tbl) + 1, ] <- c(
        "<i>Remaining variables</i>",
        "-", "-", "❌ Excluded",
        sprintf("Power limit reached (SPP drops below %d)", conf$min_samples_per_param),
        "-"
      )
      break # STOP THE LOOP!
    }

    # 3. Contrast Correlation Check (Signal Eraser)
    rho_contrast <- cor_matrices$rho[v, group_col]
    p_contrast <- cor_matrices$p[v, group_col]

    if (!is.na(rho_contrast) && !is.na(p_contrast) && abs(rho_contrast) >= conf$max_cor_rho && p_contrast <= conf$max_cor_p) {
      if (is_forced) {
        rec <- "✅ Include"; reason <- sprintf("FORCED. WARNING: Confounder with contrast (rho=%.2f)", rho_contrast)
        selected_covariates <- c(selected_covariates, v)
      } else {
        results_tbl[i, ] <- c(v, vpa_rank, dist_str, "❌ Exclude", sprintf("Confounder with contrast (rho=%.2f)", rho_contrast), "-")
        next
      }
    } else {
      # 4. Collinearity Check with Already Selected Variables
      conflict_var <- NULL
      for (sel_v in selected_covariates) {
        rho_cov <- cor_matrices$rho[v, sel_v]
        p_cov <- cor_matrices$p[v, sel_v]
        if (!is.na(rho_cov) && !is.na(p_cov) && abs(rho_cov) >= conf$max_cor_rho && p_cov <= conf$max_cor_p) {
          conflict_var <- sprintf("%s (rho=%.2f)", sel_v, rho_cov)
          break
        }
      }

      if (!is.null(conflict_var)) {
        if (is_forced) {
          rec <- "✅ Include"; reason <- sprintf("FORCED. WARNING: Redundant with %s", conflict_var)
          selected_covariates <- c(selected_covariates, v)
        } else {
          results_tbl[i, ] <- c(v, vpa_rank, dist_str, "❌ Exclude", sprintf("Redundant with included variable %s", conflict_var), "-")
          next
        }
      } else {
        # Passed all checks!
        rec <- "✅ Include"; reason <- if(is_forced) "Forced (Clean)" else "Top independent VPA hit"
        selected_covariates <- c(selected_covariates, v)
      }
    }

    current_cc_n <- temp_n
    results_tbl[i, ] <- c(v, vpa_rank, dist_str, rec, reason, "-")
  }

  results_tbl <- results_tbl %>% drop_na(Variable)

  # 5. Calculate VIF for Final Model
  final_vars <- c(group_col, selected_covariates)
  cc_df <- sub_df %>% select(Subject_ID, all_of(final_vars)) %>% drop_na()

  matrix_status <- "Full Rank (Independent)"
  if (length(selected_covariates) > 1) {
    # Generate dummy variable to trick lm() into calculating VIF of covariates
    cc_df$dummy_y <- rnorm(nrow(cc_df))
    form_str <- paste("dummy_y ~", paste(sprintf("`%s`", selected_covariates), collapse = " + "))

    vif_res <- tryCatch(car::vif(lm(as.formula(form_str), data = cc_df)), error = function(e) NULL)

    if (!is.null(vif_res)) {
      for (v in names(vif_res)) {
        # BUG FIX: car::vif returns names with backticks. Strip them so they match the table!
        v_clean <- gsub("`", "", v)

        # car::vif returns matrix for factors, we just take the first column (GVIF)
        v_val <- if(is.matrix(vif_res) || is.data.frame(vif_res)) vif_res[v, 1] else vif_res[v]
        results_tbl$VIF[results_tbl$Variable == v_clean] <- sprintf("%.2f", v_val)
      }
    } else {
      matrix_status <- "<span style='color:#c0392b; font-weight:bold;'>Rank Deficient (Aliased)</span>"
    }
  } else if (length(selected_covariates) == 1) {
    results_tbl$VIF[results_tbl$Variable == selected_covariates[1]] <- "1.00 (Single Var)"
  }

  # Return Payload (Added N_Start so wrapper can calculate percentages)
  return(list(
    status = "success",
    table = results_tbl,
    selected_vars = selected_covariates,
    meta = list(N = current_cc_n, N_Start = n_start, SPP = current_cc_n / length(final_vars), Params = length(final_vars), Matrix = matrix_status)
  ))
}

# 9.Smart Mass Correlation (Supports Continuous & Categorical) -----------------
run_smart_mass_cor <- function(data_mat, score_vec, dict_df = NULL) {
  # data_mat: variables in rows, patients in columns
  # score_vec: named vector of the continuous score

  res_list <- list()

  for (feat in rownames(data_mat)) {
    x <- data_mat[feat, names(score_vec)]

    # 1. Skip if too much missing data
    if(sum(!is.na(x)) < 10) next

    # 2. Determine Data Type
    d_type <- "Continuous" # Default for proteomics
    if (!is.null(dict_df) && feat %in% dict_df$Variable) {
      d_type <- dict_df$Class[dict_df$Variable == feat]
    }

    # 3. Route to the correct statistical test
    if (d_type == "Categorical_Binary") {
      # Wilcoxon Rank-Sum (Returns effect size r instead of Rho)
      # Note: requires numeric 0/1 encoding
      groups <- factor(x)
      if (length(levels(groups)) == 2) {
        wt <- wilcox.test(score_vec ~ groups)
        # Approximate effect size r = Z / sqrt(N)
        z <- qnorm(wt$p.value / 2)
        eff_size <- abs(z) / sqrt(length(x))
        res_list[[feat]] <- c(Feature=feat, Effect=eff_size, P_Val=wt$p.value, Method="Wilcoxon")
      }
    } else if (d_type == "Categorical_Nominal") {
      # Kruskal-Wallis for multi-level factors
      groups <- factor(x)
      if (length(levels(groups)) > 2) {
        kw <- kruskal.test(score_vec ~ groups)
        res_list[[feat]] <- c(Feature=feat, Effect=NA, P_Val=kw$p.value, Method="Kruskal-Wallis")
      }
    } else {
      # Standard Spearman for Continuous
      ct <- cor.test(x, score_vec, method="spearman", exact=FALSE)
      res_list[[feat]] <- c(Feature=feat, Effect=unname(ct$estimate), P_Val=unname(ct$p.value), Method="Spearman")
    }
  }

  # 4. Compile and FDR correct
  df <- bind_rows(lapply(res_list, function(x) as.data.frame(t(x)))) %>%
    mutate(
      Effect = as.numeric(Effect),
      P_Val = as.numeric(P_Val),
      FDR = p.adjust(P_Val, method="BH")
    ) %>%
    arrange(P_Val)

  return(df)
}

# 10. Shared Clustering Worker Engine ----------------------------------------
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

# 11. Algorithm-Aware Diagnostic & Clinical Correlation Engine -----------------
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
