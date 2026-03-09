run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  step <- "04_feature_importance_and_shap"
  log_message(context, step, "Preparing feature-importance and SHAP outputs")

  step3_dir <- file.path(context$project_root, "step3")
  legacy_models_path <- file.path(step3_dir, "models_results_v9.3.9_enhanced.csv")
  legacy_shap_candidates <- c(
    file.path(step3_dir, "plots", "SHAP_Main_OuterOOF_SVM_RFE_LDA_GlobalBar.pdf"),
    file.path(step3_dir, "plots", "SHAP_Main_OuterOOF_SVM_RFE_LDA_v9.3.9_GlobalBar.pdf")
  )

  if (file.exists(legacy_models_path) && any(file.exists(legacy_shap_candidates))) {
    log_message(context, step, "Using legacy SHAP and model summary as primary feature-importance source")

    models_df <- read.csv(legacy_models_path, check.names = FALSE)
    all_features <- unlist(strsplit(paste(models_df$Features, collapse = ";"), ";", fixed = TRUE))
    all_features <- trimws(all_features)
    all_features <- all_features[nzchar(all_features)]
    freq_table <- sort(table(all_features), decreasing = TRUE)
    freq_df <- data.frame(
      feature = names(freq_table),
      selection_frequency = as.integer(freq_table),
      stringsAsFactors = FALSE
    )

    locked_row <- models_df[models_df$Model == "SVM_RFE+LDA", , drop = FALSE]
    if (nrow(locked_row) == 0L) {
      locked_row <- models_df[1, , drop = FALSE]
    }
    locked_features <- trimws(unlist(strsplit(as.character(locked_row$Features[1]), ";", fixed = TRUE)))
    locked_features <- unique(locked_features[nzchar(locked_features)])
    locked_df <- data.frame(
      feature = locked_features,
      locked_signature_member = TRUE,
      stringsAsFactors = FALSE
    )

    importance_df <- merge(freq_df, locked_df, by = "feature", all = TRUE)
    importance_df$selection_frequency[is.na(importance_df$selection_frequency)] <- 0L
    importance_df$locked_signature_member[is.na(importance_df$locked_signature_member)] <- FALSE
    importance_df$mean_abs_shap <- NA_real_
    importance_df$permutation_auc_drop <- NA_real_
    importance_df$coefficient <- NA_real_
    importance_df <- importance_df[order(-importance_df$locked_signature_member, -importance_df$selection_frequency, importance_df$feature), , drop = FALSE]

    importance_path <- file.path(context$tables_dir, "Table_FeatureImportance.csv")
    alias_table_path <- file.path(context$tables_dir, "Table_FeatureImportance_RF.csv")
    write.csv(importance_df, importance_path, row.names = FALSE)
    write.csv(importance_df, alias_table_path, row.names = FALSE)

    legacy_shap_src <- legacy_shap_candidates[file.exists(legacy_shap_candidates)][1]
    shap_path <- file.path(context$figures_supp_dir, "SHAP_Main_OuterOOF_SVM_RFE_LDA_GlobalBar.pdf")
    file.copy(legacy_shap_src, shap_path, overwrite = TRUE)

    interp_lines <- c(
      "Feature-importance source: legacy V9.3.9 SHAP export + cross-model selection frequency summary.",
      sprintf("Locked model used for manuscript alignment: %s", as.character(locked_row$Model[1])),
      sprintf("Locked signature feature count: %d", length(locked_features)),
      "Interpretation note: SHAP bar values are preserved from the legacy artifact; machine-readable table includes stable feature membership and selection-frequency context."
    )
    write_text_file(file.path(context$tables_dir, "feature_importance_interpretation.txt"), interp_lines)

    log_message(context, step, sprintf("Feature-importance table written: %s", importance_path))
    log_message(context, step, sprintf("SHAP figure synced from legacy output: %s", shap_path))
    return(invisible(NULL))
  }

  log_message(context, step, "Legacy SHAP artifacts unavailable; falling back to canonical fastshap workflow", level = "WARN")

  obj <- readRDS(file.path(context$tables_dir, "ml_canonical_pipeline_objects.rds"))
  outer_df <- read.csv(file.path(context$tables_dir, "Nested_CV_OuterFold_Selections.csv"), check.names = FALSE)

  selected_counts <- table(unlist(strsplit(paste(outer_df$selected_features, collapse = ";"), ";", fixed = TRUE)))
  selected_counts <- selected_counts[names(selected_counts) != ""]
  selection_df <- data.frame(feature = names(selected_counts), selection_frequency = as.integer(selected_counts), stringsAsFactors = FALSE)

  coef_df <- obj$coefficients
  coef_df$abs_coefficient <- abs(coef_df$coefficient)
  coef_df <- merge(coef_df, selection_df, by = "feature", all.x = TRUE)
  coef_df$selection_frequency[is.na(coef_df$selection_frequency)] <- 0L

  baseline_prob <- predict_glmnet_prob(obj$final_fit, obj$train_x)
  baseline_auc <- roc_auc_ci(obj$train_y, baseline_prob)["auc"]
  set.seed(20260306)
  permutation_rows <- lapply(seq_len(ncol(obj$train_x)), function(i) {
    feature <- colnames(obj$train_x)[i]
    drops <- replicate(20, {
      x_perm <- obj$train_x
      x_perm[, i] <- sample(x_perm[, i])
      perm_prob <- predict_glmnet_prob(obj$final_fit, x_perm)
      baseline_auc - roc_auc_ci(obj$train_y, perm_prob)["auc"]
    })
    data.frame(feature = feature, permutation_auc_drop = mean(drops, na.rm = TRUE), stringsAsFactors = FALSE)
  })
  permutation_df <- do.call(rbind, permutation_rows)

  pred_wrapper <- function(object, newdata) {
    predict_glmnet_prob(object, as.matrix(newdata))
  }
  shap_values <- fastshap::explain(
    object = obj$final_fit,
    X = as.data.frame(obj$train_x),
    pred_wrapper = pred_wrapper,
    nsim = 64,
    adjust = TRUE
  )
  shap_importance <- data.frame(
    feature = colnames(shap_values),
    mean_abs_shap = vapply(shap_values, function(x) mean(abs(x)), numeric(1)),
    stringsAsFactors = FALSE
  )

  importance_df <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), list(coef_df, permutation_df, shap_importance))
  importance_df[is.na(importance_df)] <- 0
  importance_df <- importance_df[order(-importance_df$mean_abs_shap, -importance_df$permutation_auc_drop), , drop = FALSE]

  importance_path <- file.path(context$tables_dir, "Table_FeatureImportance.csv")
  write.csv(importance_df, importance_path, row.names = FALSE)
  alias_table_path <- file.path(context$tables_dir, "Table_FeatureImportance_RF.csv")
  write.csv(importance_df, alias_table_path, row.names = FALSE)

  plot_df <- utils::head(importance_df, 15)
  plot_df <- plot_df[!duplicated(plot_df$feature), , drop = FALSE]
  plot_df$feature <- factor(plot_df$feature, levels = rev(unique(plot_df$feature)))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = feature, y = mean_abs_shap)) +
    ggplot2::geom_col(fill = "#B85C38") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Global SHAP importance from the locked nested-CV model",
      subtitle = "Elastic-net logistic regression trained on GSE89632",
      x = NULL,
      y = "Mean absolute SHAP value"
    ) +
    ggplot2::theme_bw(base_size = 11)
  shap_path <- file.path(context$figures_supp_dir, "SHAP_Main_OuterOOF_SVM_RFE_LDA_GlobalBar.pdf")
  ggplot2::ggsave(shap_path, plot = p, width = 7.5, height = 5.5)

  if (sum(importance_df$mean_abs_shap, na.rm = TRUE) > 0) {
    dominance <- importance_df$mean_abs_shap[1] / sum(importance_df$mean_abs_shap)
  } else {
    dominance <- NA_real_
  }
  interp_lines <- c(
    sprintf("Top feature by SHAP: %s", importance_df$feature[1]),
    sprintf("Top-feature SHAP share: %.3f", dominance)
  )
  if (is.finite(dominance) && dominance > 0.80) {
    interp_lines <- c(
      interp_lines,
      "One feature dominates the model explanation. Suggested wording downgrade: describe the model as largely driven by a single marker with supporting multigene context, rather than a balanced signature."
    )
  } else {
    interp_lines <- c(interp_lines, "Feature contributions are distributed enough to support a multigene-program interpretation.")
  }
  write_text_file(file.path(context$tables_dir, "feature_importance_interpretation.txt"), interp_lines)

  log_message(context, step, sprintf("Feature-importance table written: %s", importance_path))
  log_message(context, step, sprintf("SHAP figure written: %s", shap_path))
}
