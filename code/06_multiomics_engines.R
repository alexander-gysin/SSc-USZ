# Multiomics Engines

# WRAPPERS ---------------------------------------------------------------------

wrap_investigator_correlation <- function(target_vec, target_type, olink_mat, lipid_mat, clin_mat, clin_dict) {

  # Diagnostic version of assign_direction
  assign_direction <- function(df, domain_name) {
    if(nrow(df) == 0) return(df)

    # DIAGNOSTIC INTERCEPT
    if (!"Effect" %in% colnames(df)) {
      message(sprintf("\n--- DIAGNOSTIC ALERT: %s Domain ---", domain_name))
      message("Legacy engine did not return an 'Effect' column. Available columns are:")
      message(paste(colnames(df), collapse = ", "))
      message("First 3 rows:")
      print(head(df, 3))
      message("-------------------------------------------")

      # Mock the column to prevent the loop from crashing
      df$Effect <- NA
    }

    df %>% mutate(
      Effect = as.numeric(Effect),
      Statistical_Direction = case_when(
        Method %in% c("Spearman", "Wilcoxon") & Effect > 0 ~ "Positive (Up)",
        Method %in% c("Spearman", "Wilcoxon") & Effect < 0 ~ "Negative (Down)",
        Method == "Fisher_Exact" & Effect > 1 ~ "Positive Association",
        Method == "Fisher_Exact" & Effect < 1 ~ "Negative Association",
        TRUE ~ "Unknown"
      )
    )
  }

  olink_dict <- data.frame(Variable = rownames(olink_mat), Class = "Continuous")
  res_olink <- run_smart_mass_cor(data_mat = olink_mat, score_vec = target_vec, dict_df = olink_dict)
  res_olink$Domain <- "Proteomics"

  lipid_dict <- data.frame(Variable = rownames(lipid_mat), Class = "Categorical_Binary")
  if (target_type == "Continuous") {
    res_lipid <- run_smart_mass_cor(data_mat = lipid_mat, score_vec = target_vec, dict_df = lipid_dict)
  } else {
    res_lipid <- run_fisher_cor_binary(data_mat = lipid_mat, score_vec = target_vec)
  }
  res_lipid$Domain <- "Lipidomics"

  res_clin <- run_smart_mass_cor(data_mat = clin_mat, score_vec = target_vec, dict_df = clin_dict)
  res_clin$Domain <- "Clinical"

  return(list(
    Proteomics = assign_direction(res_olink, "Proteomics"),
    Lipidomics = assign_direction(res_lipid, "Lipidomics"),
    Clinical   = assign_direction(res_clin, "Clinical")
  ))
}

wrap_highlight_plot <- function(t_vec, target_name, t_type, hl_conf, env, plot_theme, clin_dict) {
  # (Your existing wrap_highlight_plot logic remains exactly the same here)
  hl_var <- hl_conf$var
  hl_src_name <- hl_conf$source

  if(!exists(hl_src_name, envir = env)) return(list(status = "failed", msg = sprintf("Source '%s' not found.", hl_src_name)))
  src_df <- get(hl_src_name, envir = env)

  if(!hl_var %in% colnames(src_df)) return(list(status = "failed", msg = sprintf("Variable '%s' not in '%s'.", hl_var, hl_src_name)))
  if(!"Subject_ID" %in% colnames(src_df)) return(list(status = "failed", msg = sprintf("Source '%s' missing 'Subject_ID'.", hl_src_name)))

  df_target <- data.frame(Subject_ID = names(t_vec), Target_Value = as.numeric(t_vec))
  df_hl <- src_df %>% select(Subject_ID, HL_Value = all_of(hl_var))

  df_plot <- inner_join(df_target, df_hl, by = "Subject_ID") %>% drop_na()
  n_size <- nrow(df_plot)

  if(n_size < 5) return(list(status = "failed", msg = sprintf("Insufficient N after dropping NAs (N=%d).", n_size)))

  is_cont <- FALSE
  if (!is.null(clin_dict) && hl_var %in% clin_dict$Variable) {
    if (clin_dict$Class[clin_dict$Variable == hl_var] == "Continuous") is_cont <- TRUE
  }

  plot_title <- sprintf("%s vs %s (N = %d)", target_name, hl_var, n_size)

  if(is_cont) {
    df_plot$HL_Value <- as.numeric(df_plot$HL_Value)
    fit <- lm(Target_Value ~ HL_Value, data = df_plot)
    f_stat <- summary(fit)$fstatistic
    p_val <- if(!is.null(f_stat)) pf(f_stat[1], f_stat[2], f_stat[3], lower.tail=FALSE) else NA
    subtitle <- sprintf("Linear Regression: R2 = %.3f, p = %.3e", summary(fit)$r.squared, p_val)

    p <- ggplot(df_plot, aes(x = HL_Value, y = Target_Value)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", color = "darkred", se = TRUE) +
      labs(title = plot_title, subtitle = subtitle, x = hl_var, y = "Abundance") + plot_theme
  } else {
    df_plot$HL_Value <- as.factor(df_plot$HL_Value)
    if (length(levels(df_plot$HL_Value)) < 2) return(list(status = "failed", msg = sprintf("Only 1 factor level remains.")))

    if (t_type == "Continuous") {
      kw <- kruskal.test(Target_Value ~ HL_Value, data = df_plot)
      subtitle <- sprintf("Overall Test: p = %.3e", kw$p.value)
      my_comparisons <- combn(levels(df_plot$HL_Value), 2, simplify = FALSE)

      p <- ggplot(df_plot, aes(x = HL_Value, y = Target_Value, fill = HL_Value)) +
        geom_violin(alpha = 0.7, width = 0.8) +
        geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
        scale_fill_manual(values = get_project_colors(levels(df_plot$HL_Value))) +
        ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.01) +
        labs(title = plot_title, subtitle = subtitle, x = hl_var, y = "Abundance") + plot_theme + theme(legend.position = "none")
    } else {
      p <- ggplot(df_plot, aes(x = HL_Value, fill = factor(Target_Value))) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("0" = "grey80", "1" = COLOR_HIGHLIGHT)) +
        labs(title = plot_title, subtitle = "Binary Distribution Profile", x = hl_var, y = "Proportion", fill = "Detected") + plot_theme
    }
  }
  return(list(status = "success", plot = p, p_name = sprintf("%s_vs_%s", target_name, hl_var)))
}


# WORKERS -----------------------------------------------------


run_fisher_cor_binary <- function(data_mat, score_vec) {
  # (Your existing run_fisher_cor_binary remains exactly the same here)
  res_list <- list()
  for (feat in rownames(data_mat)) {
    x <- data_mat[feat, names(score_vec)]
    if(sum(!is.na(x)) < 10) next
    tbl <- table(Factor_A = factor(score_vec, levels=c(0,1)), Factor_B = factor(x, levels=c(0,1)))
    if (all(dim(tbl) == c(2, 2))) {
      ft <- fisher.test(tbl)
      eff_size <- unname(ft$estimate)
      p_val <- ft$p.value
    } else {
      eff_size <- NA
      p_val <- NA
    }
    res_list[[feat]] <- c(Feature=feat, Effect=eff_size, P_Val=p_val, Method="Fisher_Exact")
  }
  df <- bind_rows(lapply(res_list, function(x) as.data.frame(t(x)))) %>%
    mutate(Effect = as.numeric(Effect), P_Val = as.numeric(P_Val), FDR = p.adjust(P_Val, method="BH")) %>%
    filter(!is.na(P_Val)) %>% arrange(P_Val)
  return(df)
}
