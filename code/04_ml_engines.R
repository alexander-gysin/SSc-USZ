# ==============================================================================
# TIER 1: ORCHESTRATORS (Wrappers)
# ==============================================================================

#' Wrapper: ML Data Preparation
#' Routes omics types, applies hybrid sparsity logic for lipidomics,
#' and standardizes matrices (Z-scoring).
wrap_ml_data_prep <- function(mat, clin_df, config) {

  # 1. Align Clinical and Omics Data
  clin_sub <- clin_df %>%
    dplyr::filter(.data[[config$test_var]] %in% c(config$target_groups, config$ref_groups)) %>%
    dplyr::filter(!is.na(.data[[config$test_var]]))

  common_samples <- intersect(colnames(mat), clin_sub$Subject_ID)

  if (length(common_samples) < 10) {
    stop("Insufficient overlapping samples (< 10) found between clinical and omics data.")
  }

  mat_sub <- mat[, common_samples, drop = FALSE]
  clin_sub <- clin_sub[match(common_samples, clin_sub$Subject_ID), ]

  # 2. Define Y (Target vs Reference)
  Y <- ifelse(clin_sub[[config$test_var]] %in% config$target_groups, "Target", "Reference")
  Y <- factor(Y, levels = c("Reference", "Target")) # Ensure specific factor leveling for pROC

  # 3. Define X (Features) based on Omics Topology
  X_raw <- t(mat_sub)

  if (config$omics_type == "lipidomics") {

    # Hybrid Matrix Logic: 0/1 for sparse, scaled continuous for dense
    missing_rates <- colMeans(is.na(X_raw) | X_raw == 0)

    sparse_cols <- names(missing_rates[missing_rates >= 0.30])
    dense_cols <- names(missing_rates[missing_rates < 0.30])

    X_hybrid <- matrix(nrow = nrow(X_raw), ncol = ncol(X_raw), dimnames = dimnames(X_raw))

    # Binarize sparse features
    if (length(sparse_cols) > 0) {
      X_hybrid[, sparse_cols] <- ifelse(is.na(X_raw[, sparse_cols]) | X_raw[, sparse_cols] == 0, 0, 1)
    }

    # Scale dense features (Z-Score)
    if (length(dense_cols) > 0) {
      X_hybrid[, dense_cols] <- scale(X_raw[, dense_cols, drop = FALSE])
    }

    X <- X_hybrid

  } else if (config$omics_type == "proteomics") {

    # Pure continuous logic: Scale to Z-scores
    # Note: Assuming 100% complete as discussed. If imputation is needed, it goes here.
    X <- scale(X_raw)
  }

  return(list(X = as.matrix(X), Y = Y))
}


#' Wrapper: Universal ML Engine
#' Handles Stratified Nested CV, applies SMOTE on training folds,
#' dynamically calls requested algorithm, and routes outputs to evaluation workers.
wrap_ml_engine <- function(X, Y, config, global_config, algorithm, job_name, out_dir, assets_dir, project_colors_func) {

  # 1. Set up Stratified Outer Folds
  set.seed(42)
  outer_folds <- caret::createFolds(Y, k = global_config$cv$outer_folds, returnTrain = TRUE)

  fold_results <- list()
  feature_weights <- list()

  # 2. Nested CV Loop
  for (i in seq_along(outer_folds)) {
    train_idx <- outer_folds[[i]]

    X_train <- X[train_idx, , drop = FALSE]
    Y_train <- Y[train_idx]
    X_test <- X[-train_idx, , drop = FALSE]
    Y_test <- Y[-train_idx]

    # 3. Apply SMOTE strictly to training data to prevent leakage
    if (isTRUE(config$use_smote) && config$omics_type == "proteomics") {
      # SMOTE requires numeric class (0/1) mapping
      y_train_num <- ifelse(Y_train == "Target", 1, 0)
      smote_df <- data.frame(X_train, class = y_train_num)

      smote_res <- smotefamily::SMOTE(smote_df[, -ncol(smote_df)], smote_df$class)

      X_train <- as.matrix(smote_res$data[, -ncol(smote_res$data)])
      Y_train <- factor(ifelse(smote_res$data$class == 1, "Target", "Reference"), levels = c("Reference", "Target"))
    }

    # 4. Route to Specific Algorithm Worker
    hyperparams <- config$hyperparameters[[algorithm]]

    if (algorithm == "lasso") {
      res <- worker_fit_lasso(X_train, Y_train, X_test, hyperparams, global_config$cv$inner_folds)
    } else if (algorithm == "rf") {
      res <- worker_fit_rf(X_train, Y_train, X_test, hyperparams)
    } else if (algorithm == "gb") {
      res <- worker_fit_gb(X_train, Y_train, X_test, hyperparams, global_config$cv$inner_folds)
    } else if (algorithm == "plsda") {
      res <- worker_fit_plsda(X_train, Y_train, X_test, hyperparams)
    }

    # 5. Store Fold Predictions and Weights
    fold_results[[i]] <- data.frame(Subject_ID = rownames(X_test), Fold = i, True_Y = Y_test, Pred_Prob = res$preds)

    if (nrow(res$weights) > 0) {
      feature_weights[[i]] <- res$weights %>% dplyr::mutate(Fold = i)
    }
  }

  # 6. Aggregate Results across all folds
  all_preds <- do.call(rbind, fold_results)
  all_weights <- do.call(rbind, feature_weights)

  # 7. Call Evaluation Workers
  metrics_df <- worker_calc_metrics(all_preds$True_Y, all_preds$Pred_Prob, algorithm)
  p_roc <- worker_plot_roc(all_preds, algorithm, job_name, project_colors_func)

  # Calculate Feature Stability
  if (!is.null(all_weights) && nrow(all_weights) > 0) {
    stability_df <- all_weights %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarize(
        Selection_Frequency = sum(Importance > 0) / global_config$cv$outer_folds,
        Mean_Importance = mean(Importance),
        .groups = "drop"
      ) %>%
      dplyr::filter(Selection_Frequency >= global_config$feature_stability_threshold) %>%
      dplyr::arrange(desc(Mean_Importance))

    p_stab <- worker_plot_stability(stability_df, algorithm, job_name, project_colors_func)
  } else {
    stability_df <- data.frame(Feature = character(), Selection_Frequency = numeric(), Mean_Importance = numeric())
    p_stab <- NULL
  }

  # 8. Export Artifacts
  write.csv(metrics_df, file.path(out_dir, paste0(job_name, "_", algorithm, "_metrics.csv")), row.names = FALSE)
  write.csv(stability_df, file.path(out_dir, paste0(job_name, "_", algorithm, "_stable_features.csv")), row.names = FALSE)

  ggsave(file.path(out_dir, paste0(job_name, "_", algorithm, "_ROC.png")), plot = p_roc, width = 6, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(assets_dir, paste0(job_name, "_", algorithm, "_ROC.png")), plot = p_roc, width = 6, height = 6, dpi = 300, bg = "white")

  if (!is.null(p_stab)) {
    ggsave(file.path(out_dir, paste0(job_name, "_", algorithm, "_Stability.png")), plot = p_stab, width = 8, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(assets_dir, paste0(job_name, "_", algorithm, "_Stability.png")), plot = p_stab, width = 8, height = 6, dpi = 300, bg = "white")
  }

  return(list(metrics = metrics_df, roc_plot = p_roc, stability_plot = p_stab, stability_df = stability_df))
}


# ==============================================================================
# TIER 2: EXECUTORS (Model Workers)
# ==============================================================================

#' Worker: Fit LASSO / Elastic Net
#' Inner CV handles Lambda tuning. Returns predictions and extracted weights.
worker_fit_lasso <- function(X_train, Y_train, X_test, params, inner_folds) {

  y_num <- ifelse(Y_train == "Target", 1, 0)

  # Calculate Class Weights to penalize minority class errors
  wts <- ifelse(y_num == 1, 1 / sum(y_num == 1), 1 / sum(y_num == 0))
  wts <- wts / mean(wts) # Normalize sum to N

  # Fit with inner CV for lambda
  cv_fit <- glmnet::cv.glmnet(
    x = X_train, y = y_num,
    family = "binomial",
    alpha = params$alpha,
    weights = wts,
    nfolds = inner_folds,
    type.measure = "auc"
  )

  # Predict on blinded test fold
  preds <- predict(cv_fit, newx = X_test, s = "lambda.min", type = "response")

  # Extract feature coefficients
  coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
  weights <- data.frame(Feature = rownames(coefs)[-1], Importance = abs(coefs[-1, 1])) %>%
    dplyr::filter(Importance > 0)

  return(list(preds = as.numeric(preds), weights = weights))
}


#' Worker: Fit Random Forest
#' Handles class imbalance natively via sampling parameters.
worker_fit_rf <- function(X_train, Y_train, X_test, params) {

  # Downsample majority class strictly within training boundaries for balance
  min_class_size <- min(table(Y_train))

  rf_fit <- randomForest::randomForest(
    x = X_train, y = Y_train,
    ntree = params$ntree,
    strata = Y_train,
    sampsize = c(min_class_size, min_class_size),
    importance = TRUE
  )

  preds <- predict(rf_fit, newdata = X_test, type = "prob")[, "Target"]

  imp <- randomForest::importance(rf_fit)
  weights <- data.frame(Feature = rownames(imp), Importance = imp[, "MeanDecreaseGini"]) %>%
    dplyr::filter(Importance > 0)

  return(list(preds = as.numeric(preds), weights = weights))
}


#' Worker: Fit Gradient Boosting (XGBoost)
#' Includes Early Stopping mapped to inner folds to prevent severe overfitting.
worker_fit_gb <- function(X_train, Y_train, X_test, params, inner_folds) {

  y_train_num <- ifelse(Y_train == "Target", 1, 0)
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train_num)
  dtest <- xgboost::xgb.DMatrix(data = X_test)

  # Scale pos weight for class imbalance
  spw <- sum(y_train_num == 0) / sum(y_train_num == 1)

  xgb_params <- list(
    objective = "binary:logistic",
    max_depth = params$max_depth,
    eta = params$eta,
    scale_pos_weight = spw,
    eval_metric = "auc"
  )

  # Inner CV purely for Early Stopping
  cv_xgb <- xgboost::xgb.cv(
    params = xgb_params,
    data = dtrain,
    nrounds = params$nrounds,
    nfold = inner_folds,
    showsd = FALSE,
    stratified = TRUE,
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose = 0
  )

  best_iter <- cv_xgb$best_iteration

  # Train final model on entire Train fold using best iteration
  fit <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = best_iter, verbose = 0)

  preds <- predict(fit, dtest)
  imp <- xgboost::xgb.importance(model = fit)

  # Catch edge case where GB learns absolutely nothing
  if (nrow(imp) > 0) {
    weights <- data.frame(Feature = imp$Feature, Importance = imp$Gain)
  } else {
    weights <- data.frame(Feature = character(), Importance = numeric())
  }

  return(list(preds = as.numeric(preds), weights = weights))
}


#' Worker: Fit PLS-DA
#' Extracts optimal components and approximates VIP.
worker_fit_plsda <- function(X_train, Y_train, X_test, params) {

  # Suppress caret verbosity during training
  fit <- suppressWarnings(
    caret::plsda(X_train, Y_train, ncomp = params$ncomp)
  )

  preds <- predict(fit, X_test, type = "prob")[, "Target", 1]

  # VIP approximation via caret varImp
  imp <- caret::varImp(fit)$importance
  weights <- data.frame(Feature = rownames(imp), Importance = imp[, "Target"]) %>%
    dplyr::filter(Importance > 0)

  return(list(preds = as.numeric(preds), weights = weights))
}


# ==============================================================================
# TIER 3: EVALUATORS (Utility Workers)
# ==============================================================================

#' Worker: Calculate Unified Metrics
#' Computes AUROC and Youden's Index for optimal specificty/sensitivity
worker_calc_metrics <- function(true_y, pred_prob, alg_name) {

  roc_obj <- pROC::roc(response = true_y, predictor = pred_prob, levels = c("Reference", "Target"), direction = "<", quiet = TRUE)
  auroc <- as.numeric(pROC::auc(roc_obj))

  # Optimal threshold by Youden's J statistic
  coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"), best.method = "youden", transpose = FALSE)

  # Guard against multiple thresholds returned
  df_metrics <- data.frame(
    Algorithm = toupper(alg_name),
    AUROC = round(auroc, 3),
    Sensitivity = round(coords$sensitivity[1], 3),
    Specificity = round(coords$specificity[1], 3),
    Optimal_Threshold = round(coords$threshold[1], 3)
  )

  return(df_metrics)
}


#' Worker: Plot Standardized ROC
worker_plot_roc <- function(pred_df, alg_name, job_name, project_colors_func) {

  roc_obj <- pROC::roc(response = pred_df$True_Y, predictor = pred_df$Pred_Prob, levels = c("Reference", "Target"), direction = "<", quiet = TRUE)
  auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)

  roc_data <- data.frame(TPR = roc_obj$sensitivities, FPR = 1 - roc_obj$specificities)

  # Extract dynamic target color (assumes "Target" maps roughly to your active config)
  line_color <- project_colors_func("Active")

  p <- ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = line_color, linewidth = 1.2) +
    ggplot2::geom_abline(linetype = "dashed", color = "grey50") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = sprintf("ROC: %s", toupper(alg_name)),
      subtitle = sprintf("AUROC = %s", auc_val),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  return(p)
}


#' Worker: Plot Feature Stability
worker_plot_stability <- function(stability_df, alg_name, job_name, project_colors_func) {

  if (nrow(stability_df) == 0) return(NULL)

  top_features <- head(stability_df, 20)
  bar_color <- project_colors_func("Non-Active")

  p <- ggplot2::ggplot(top_features, ggplot2::aes(x = reorder(Feature, Mean_Importance), y = Mean_Importance, fill = Selection_Frequency)) +
    ggplot2::geom_bar(stat = "identity", color = "grey20", size = 0.2) +
    ggplot2::scale_fill_gradient(low = "grey90", high = bar_color, limits = c(0, 1), name = "CV Stability") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = sprintf("Stable Signatures: %s", toupper(alg_name)),
      x = NULL,
      y = "Mean CV Importance"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  return(p)
}
