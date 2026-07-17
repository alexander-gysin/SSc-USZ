# Multiomics Engines

# WRAPPERS ---------------------------------------------------------------------

wrap_investigator_correlation <- function(target_vec, target_type, olink_mat, lipid_mat, clin_mat, clin_dict) {
  #
  # target_vec: named numeric vector
  # target_type: "Continuous" (Proteomics) or "Categorical_Binary" (Lipidomics)

  results <- list()

  # 1. vs Proteomics (Continuous Matrix)
  olink_dict <- data.frame(Variable = rownames(olink_mat), Class = "Continuous")

  # Legacy engine will naturally use Spearman if target is continuous OR binary
  res_olink <- run_smart_mass_cor(data_mat = olink_mat, score_vec = target_vec, dict_df = olink_dict)
  res_olink$Domain <- "Proteomics"
  results$Proteomics <- res_olink

  # 2. vs Lipidomics (Binary Matrix)
  lipid_dict <- data.frame(Variable = rownames(lipid_mat), Class = "Categorical_Binary")

  if (target_type == "Continuous") {
    # Continuous Target vs Binary Matrix -> Legacy Wilcoxon
    res_lipid <- run_smart_mass_cor(data_mat = lipid_mat, score_vec = target_vec, dict_df = lipid_dict)
  } else {
    # Binary Target vs Binary Matrix -> New Fisher Worker
    res_lipid <- run_fisher_cor_binary(data_mat = lipid_mat, score_vec = target_vec)
  }
  res_lipid$Domain <- "Lipidomics"
  results$Lipidomics <- res_lipid

  # 3. vs Clinical (Mixed Matrix)
  # Legacy engine uses clin_dict to route Spearman vs Wilcoxon vs Kruskal
  res_clin <- run_smart_mass_cor(data_mat = clin_mat, score_vec = target_vec, dict_df = clin_dict)
  res_clin$Domain <- "Clinical"
  results$Clinical <- res_clin

  return(results)
}

# WORKERS -----------------------------------------------------

run_fisher_cor_binary <- function(data_mat, score_vec) {
  # data_mat: binary matrix (variables in rows, patients in columns)
  # score_vec: named binary vector of the target
  res_list <- list()

  for (feat in rownames(data_mat)) {
    x <- data_mat[feat, names(score_vec)]

    # 1. Skip if too much missing data (though binary should be complete, safety first)
    if(sum(!is.na(x)) < 10) next

    # 2. Fisher's Exact Test
    tbl <- table(Factor_A = factor(score_vec, levels=c(0,1)),
                 Factor_B = factor(x, levels=c(0,1)))

    # Only run if matrix has 2x2 dimensions (some features might be all 0s or all 1s in the subset)
    if (all(dim(tbl) == c(2, 2))) {
      ft <- fisher.test(tbl)
      # Extract Odds Ratio as effect size
      eff_size <- unname(ft$estimate)
      p_val <- ft$p.value
    } else {
      eff_size <- NA
      p_val <- NA
    }

    res_list[[feat]] <- c(Feature=feat, Effect=eff_size, P_Val=p_val, Method="Fisher_Exact")
  }

  # 3. Compile and FDR correct
  df <- bind_rows(lapply(res_list, function(x) as.data.frame(t(x)))) %>%
    mutate(
      Effect = as.numeric(Effect),
      P_Val = as.numeric(P_Val),
      FDR = p.adjust(P_Val, method="BH")
    ) %>%
    filter(!is.na(P_Val)) %>%
    arrange(P_Val)

  return(df)
}
