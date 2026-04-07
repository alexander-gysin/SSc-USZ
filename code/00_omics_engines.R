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
# 1. PCA ENGINE
# ------------------------------------------------------------------------------
run_pca_engine <- function(omics_mat, clin_df, pca_conf, project_colors) {
  if (!is.null(pca_conf$subset_to_contrast)) {
    target_col <- pca_conf$subset_to_contrast$group_col
    target_lvls <- pca_conf$subset_to_contrast$target_groups
    clin_df <- clin_df %>% filter(!!sym(target_col) %in% target_lvls)
  }
  valid_barcodes <- intersect(colnames(omics_mat), clin_df$Subject_ID)
  clin_df <- clin_df %>% filter(Subject_ID %in% valid_barcodes)
  mat_sub <- omics_mat[, clin_df$Subject_ID, drop = FALSE]

  pca_res <- FactoMineR::PCA(t(mat_sub), scale.unit = TRUE, graph = FALSE)
  generated_plots <- list()

  generated_plots[["scree"]] <- fviz_eig(pca_res, addlabels = TRUE) + theme_project_base() + labs(title = sprintf("%s: Variance Explained (Scree)", pca_conf$title))
  generated_plots[["loadings"]] <- fviz_pca_var(pca_res, select.var = list(contrib = 30), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) + theme_project_base() + labs(title = sprintf("%s: PCA Loadings (Top 30)", pca_conf$title))
  generated_plots[["pc1_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 1, top = 10, fill = "grey30", color = "grey30") + theme_project_base() + coord_flip() + labs(title = "Top PC1 Drivers")
  generated_plots[["pc2_drivers"]] <- fviz_contrib(pca_res, choice = "var", axes = 2, top = 10, fill = "grey30", color = "grey30") + theme_project_base() + coord_flip() + labs(title = "Top PC2 Drivers")

  for (c_var in pca_conf$color_vars) {
    if (!(c_var %in% colnames(clin_df))) next
    dynamic_title <- sprintf("%s (Colored by %s)", pca_conf$title, c_var)
    if (is.numeric(clin_df[[c_var]])) {
      p <- fviz_pca_ind(pca_res, label = "none", col.ind = clin_df[[c_var]], legend.title = c_var) + scale_color_viridis_c(option = "viridis") + theme_project_base() + labs(title = dynamic_title)
    } else {
      clin_df[[c_var]] <- as.factor(clin_df[[c_var]])
      p <- fviz_pca_ind(pca_res, label = "none", habillage = clin_df[[c_var]], palette = get_project_colors(levels(clin_df[[c_var]])), addEllipses = TRUE) + theme_project_base() + labs(title = dynamic_title)
    }
    generated_plots[[paste0("overlay_", c_var)]] <- p
  }
  return(generated_plots)
}

# ------------------------------------------------------------------------------
# 2. LIMMA & VOLCANO ENGINE (Individual Cutoffs)
# ------------------------------------------------------------------------------
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
    formula_str <- paste("~ 0 + limma_group +", paste(dea_conf$covariates, collapse = " + "))
  } else {
    formula_str <- "~ 0 + limma_group"
  }

  design <- model.matrix(as.formula(formula_str), data = clin_sub)
  colnames(design)[1:2] <- c("Reference", "Target")

  fit <- lmFit(mat_sub, design)
  contrast_mat <- makeContrasts(Target - Reference, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_mat))

  # INDIVIDUAL CUTOFFS APPLIED HERE
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
    scale_color_manual(values = c("Down-regulated" = ref_color, "Not Significant" = "grey85", "Up-regulated" = target_color)) +
    geom_text_repel(data = top_labels, aes(label = Feature), size = 3, color = "black", box.padding = 0.5, max.overlaps = Inf) +
    geom_vline(xintercept = c(-dea_conf$dea_fc_cutoff, dea_conf$dea_fc_cutoff), linetype = "dotted") +
    geom_hline(yintercept = -log10(dea_conf$dea_p_cutoff), linetype = "dashed", color = "darkred") +
    theme_project_base() +
    labs(title = sprintf("Volcano: %s", dea_conf$title), x = "Log2 Fold Change", y = "-log10(FDR)")

  return(list(data = dea_res, volcano = p_volc))
}

# ------------------------------------------------------------------------------
# 3. GSEA ENGINE (Fixed Title Cutoff)
# ------------------------------------------------------------------------------
run_gsea_engine <- function(dea_res, go_db, omics_config, title) {
  ranked_vec <- setNames(dea_res$t, dea_res$Feature) %>% sort(decreasing = TRUE)
  set.seed(42)
  gsea_res <- GSEA(geneList = ranked_vec, TERM2GENE = go_db, pvalueCutoff = 1, minGSSize = omics_config$gsea_min_size, verbose = FALSE)
  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(NULL)

  gsea_res@result$Description <- gsub("GOBP_", "", gsea_res@result$Description) %>% gsub("_", " ", .)
  res_df <- as.data.frame(gsea_res) %>% mutate(Status = ifelse(NES > 0, "Up-regulated", "Down-regulated"))

  warning_tag <- if(sum(res_df$p.adjust < omics_config$gsea_p_cutoff) == 0) "\n(Exploratory: No paths passed FDR)" else ""
  plot_df <- res_df %>% group_by(Status) %>% slice_max(abs(NES), n = 10) %>% ungroup()

  # FIXED TITLE: Wrap long contrast titles so they physically cannot overrun the right margin
  main_title <- str_wrap(sprintf("GSEA: %s", title), width = 60)
  full_title <- paste0(main_title, warning_tag)

  p_dot <- ggplot(plot_df, aes(x = NES, y = reorder(Description, NES), color = p.adjust, size = setSize)) +
    geom_point() +
    scale_color_gradient(low = "firebrick3", high = "navy") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme_project_base() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0)
    ) +
    labs(title = full_title, y = NULL)

  # Apply the same wrapping logic to the ridgeplot
  ridge_title <- str_wrap(sprintf("Pathways: %s", title), width = 60)
  p_ridge <- if(nrow(res_df) >= 2) ridgeplot(pairwise_termsim(gsea_res), showCategory = 15) + theme_project_base() + labs(title = ridge_title) else NULL

  return(list(data = res_df, dotplot = p_dot, ridgeplot = p_ridge))
}

# ------------------------------------------------------------------------------
# 4. CROSS-CONTRAST ENGINE (Inherits Significance)
# ------------------------------------------------------------------------------
run_cross_contrast_engine <- function(dea_x, dea_y, conf, x_label, y_label, color_x, color_y, color_shared) {

  x_sig_name <- sprintf("%s Specific", x_label)
  y_sig_name <- sprintf("%s Specific", y_label)

  temporal_shifts <- full_join(
    dea_x %>% select(Feature, logFC_X = logFC, FDR_X = adj.P.Val, Sig_Status_X = Significance),
    dea_y %>% select(Feature, logFC_Y = logFC, FDR_Y = adj.P.Val, Sig_Status_Y = Significance),
    by = "Feature"
  ) %>%
    mutate(
      # Safely inherit the pre-calculated significance from the Limma results
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

  labels_shared <- sig_points %>%
    filter(Signature == "Shared Drivers") %>%
    group_by(Quadrant) %>%
    slice_max(order_by = label_score, n = top_n_shar, with_ties = FALSE) %>%
    ungroup()

  labels_specific <- sig_points %>%
    filter(Signature %in% c(x_sig_name, y_sig_name)) %>%
    group_by(Signature) %>%
    slice_max(order_by = label_score, n = top_n_spec, with_ties = FALSE) %>%
    ungroup()

  top_labels <- bind_rows(labels_shared, labels_specific)

  sig_points <- sig_points %>% mutate(repel_label = ifelse(Feature %in% top_labels$Feature, Feature, ""))

  method_caption <- sprintf("Labeling: Top %d Shared Drivers (per quadrant) + Top %d Specific Drivers (Ranked by Pi-Score)", top_n_shar, top_n_spec)
  color_mapping <- setNames(c(color_shared, color_x, color_y), c("Shared Drivers", x_sig_name, y_sig_name))
  clean_title <- if(grepl("^Trajectory", conf$title)) paste("DEA", conf$title) else sprintf("DEA Trajectory: %s", conf$title)

  p_shift <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = 0, color = "grey40") + geom_vline(xintercept = 0, color = "grey40") +
    geom_point(data = sig_points, aes(x = logFC_X, y = logFC_Y, color = Signature), alpha = 0.8, size = 2.5) +
    geom_text_repel(
      data = sig_points, aes(x = logFC_X, y = logFC_Y, label = repel_label),
      color = "black", size = 3.5, box.padding = 0.8, point.padding = 0, min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_color_manual(values = color_mapping) +
    theme_project_base() +
    theme(plot.caption = element_text(hjust = 0, face = "italic", color = "grey40", size = 9, margin = margin(t = 15))) +
    labs(title = clean_title, x = paste(x_label, "(Log2 Fold Change)"), y = paste(y_label, "(Log2 Fold Change)"), caption = method_caption) +
    coord_fixed(ratio = 1)

  return(p_shift)
}

# ------------------------------------------------------------------------------
# 5. HEATMAP ENGINE
# ------------------------------------------------------------------------------
run_heatmap_engine <- function(mat, clin, dea_res, title, n_top = 50, split_by_group = FALSE) {
  valid_ids <- intersect(clin$Subject_ID, colnames(mat))
  clin_sub <- clin %>% filter(Subject_ID %in% valid_ids)

  top_feats <- dea_res %>% mutate(pi = abs(logFC) * -log10(adj.P.Val)) %>% slice_max(pi, n = n_top) %>% pull(Feature)
  plot_mat <- t(scale(t(mat[top_feats, clin_sub$Subject_ID])))

  grp_cols <- get_project_colors(as.character(unique(clin_sub$cohort_group)))
  sex_cols <- setNames(c(COLOR_FEMALE, COLOR_MALE), c("Female", "Male"))
  age_fun  <- colorRamp2(c(min(clin_sub$Age, na.rm=T), max(clin_sub$Age, na.rm=T)), c(COLOR_AGE_LOW, COLOR_AGE_HIGH))

  ha <- HeatmapAnnotation(Group = clin_sub$cohort_group, Sex = clin_sub$Sex, Age = clin_sub$Age,
                          col = list(Group = grp_cols, Sex = sex_cols, Age = age_fun))

  Heatmap(plot_mat, name = "Z-Score", column_title = sprintf("%s: Top %d Proteins", title, n_top),
          top_annotation = ha, show_column_names = FALSE, column_split = if(split_by_group) clin_sub$cohort_group else NULL,
          col = colorRamp2(c(-2, 0, 2), c(HM_Z_LOW, HM_Z_MID, HM_Z_HIGH)))
}

# ------------------------------------------------------------------------------
# 6. LASSO ENGINE (Comprehensive Naming & Synced Colors)
# ------------------------------------------------------------------------------
run_lasso_engine <- function(mat, clin, target_col, target_groups, title) {
  clin_sub <- clin %>% filter(!!sym(target_col) %in% target_groups)
  x <- t(mat[, clin_sub$Subject_ID])
  y <- factor(clin_sub[[target_col]])

  set.seed(42)
  cv_fit <- cv.glmnet(x, y, family = "binomial", type.measure = "auc")

  tmp_coeffs <- coef(cv_fit, s = "lambda.min")
  df_coeffs <- data.frame(Feature = rownames(tmp_coeffs), Beta = as.numeric(tmp_coeffs)) %>%
    filter(Beta != 0, Feature != "(Intercept)") %>% arrange(desc(abs(Beta)))

  # Dynamic Coloring: Beta > 0 (Up/Target) gets HM_Z_HIGH (Red), Beta < 0 (Down/Ref) gets HM_Z_LOW (Blue)
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

# ------------------------------------------------------------------------------
# 8. GSEA TRAJECTORY ENGINE (Fixed Title Cutoff)
# ------------------------------------------------------------------------------
run_pathway_trajectory_engine <- function(gsea_x, gsea_y, conf, x_lab, y_lab) {
  merged <- inner_join(
    gsea_x %>% select(Description, NES_X = NES),
    gsea_y %>% select(Description, NES_Y = NES),
    by = "Description"
  ) %>% mutate(magnitude = abs(NES_X) + abs(NES_Y))

  n_limit <- if (!is.null(conf$top_n_pathways)) conf$top_n_pathways else 20

  ggplot(merged, aes(x = NES_X, y = NES_Y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, color = "grey40") + geom_hline(yintercept = 0, color = "grey40") +
    geom_point(aes(color = magnitude), alpha = 0.7, size = 3) +
    geom_text_repel(data = slice_max(merged, magnitude, n = n_limit), aes(label = Description), size = 3) +
    scale_color_viridis_c(option = "magma", end = 0.8) +
    theme_project_base() +
    theme(plot.title.position = "plot") + # Ensures title alignment to the full canvas
    labs(title = str_wrap(sprintf("Pathway Trajectory: %s", conf$title), 50),
         x = paste(x_lab, "(NES)"), y = paste(y_lab, "(NES)"),
         color = "Shift Magnitude")
}

# ------------------------------------------------------------------------------
# 9. VENN DIAGRAM ENGINE (Inherits Significance)
# ------------------------------------------------------------------------------
run_venn_engine <- function(dea_x, dea_y, name_x, name_y, title, color_x, color_y, color_shared) {

  # Directly pull features that were already flagged as significant
  sig_x <- dea_x %>% filter(Significance != "Not Significant") %>% pull(Feature)
  sig_y <- dea_y %>% filter(Significance != "Not Significant") %>% pull(Feature)

  only_x <- length(setdiff(sig_x, sig_y))
  only_y <- length(setdiff(sig_y, sig_x))
  shared <- length(intersect(sig_x, sig_y))

  if((only_x + only_y + shared) == 0) return(NULL)

  theta <- seq(0, 2 * pi, length.out = 100)
  circle_x <- data.frame(x = -0.5 + cos(theta), y = sin(theta))
  circle_y <- data.frame(x = 0.5 + cos(theta), y = sin(theta))

  p_venn <- ggplot() +
    geom_polygon(data = circle_x, aes(x, y), fill = color_x, alpha = 0.3, color = color_x, linewidth = 1.2) +
    geom_polygon(data = circle_y, aes(x, y), fill = color_y, alpha = 0.3, color = color_y, linewidth = 1.2) +
    annotate("text", x = -0.5, y = 1.3, label = name_x, size = 5, fontface = "bold", color = color_x) +
    annotate("text", x = 0.5, y = -1.3, label = name_y, size = 5, fontface = "bold", color = color_y) +
    annotate("text", x = -0.9, y = 0, label = only_x, size = 6, fontface = "bold") +
    annotate("text", x = 0.9, y = 0, label = only_y, size = 6, fontface = "bold") +
    annotate("text", x = 0, y = 0, label = shared, size = 7, fontface = "bold", color = color_shared) +
    theme_void() +
    coord_fixed(clip = "off", xlim = c(-1.8, 1.8), ylim = c(-1.5, 1.5)) +
    labs(title = sprintf("DEA Venn: %s", title)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5), margin = margin(b = 15)),
      plot.margin = margin(10, 10, 10, 10)
    )

  return(p_venn)
}
