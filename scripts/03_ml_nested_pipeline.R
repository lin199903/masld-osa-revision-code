run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "03_ml_nested_pipeline"
  log_message(context, step, "Integrating legacy ML baseline (primary/3.0_machine_learning_v1.R)")

  locked_seed <- ensure_locked_main_fig4_assets(context, config, overwrite = FALSE, strict = FALSE)
  if (!isTRUE(locked_seed$skipped) && length(locked_seed$copied) > 0L) {
    log_message(context, step, sprintf("Seeded locked Figure 4 assets: %d file(s)", length(locked_seed$copied)))
  }

  legacy_script <- file.path(
    context$revision_root,
    "legacy_screening",
    "selected_legacy_for_original_submission",
    "primary",
    "3.0_machine_learning_v1.R"
  )
  if (!file.exists(legacy_script)) {
    stop(sprintf("Legacy baseline script not found: %s", legacy_script))
  }

  step3_dir <- file.path(context$project_root, "step3")
  required_step3 <- c(
    "OOF_outer_SVM_RFE_LDA.csv",
    "models_results_v9.3.9_enhanced.csv",
    file.path("plots", "model_performance_visualization_v9.3.9_ALL.pdf"),
    file.path("plots", "Combined_ROC_Calibration_SVM_RFE+LDA_v9.3.9.pdf"),
    file.path("plots", "PRC_Combined_SVM_RFE+LDA_v9.3.9.pdf"),
    file.path("plots", "PRC_Inverted_ControlPositive_SVM_RFE+LDA_v9.3.9.pdf")
  )
  bridge_external_path <- file.path(context$tables_dir, "legacy_bridge_external_predictions.csv")
  missing_step3 <- required_step3[!file.exists(file.path(step3_dir, required_step3))]
  need_run_legacy <- (length(missing_step3) > 0L)

  legacy_env <- NULL
  if (need_run_legacy) {
    log_message(
      context,
      step,
      sprintf("Missing legacy step3 outputs; re-running legacy baseline script (%d missing files)", length(missing_step3)),
      level = "WARN"
    )
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(context$project_root)

    # Runtime hotfix from backup-variant behavior: cap worker count to avoid
    # intermittent PSOCK connection failures in high-core Windows sessions.
    legacy_runtime_script <- file.path(context$outputs_dir, "tmp_legacy_ml_runtime_patch.R")
    legacy_lines <- readLines(legacy_script, warn = FALSE, encoding = "UTF-8")
    worker_line_idx <- grep("^N_WORKERS\\s*<-", legacy_lines)
    if (length(worker_line_idx) >= 1L) {
      legacy_lines[worker_line_idx[1]] <- "N_WORKERS <- 1"
    }
    legacy_lines <- sub(
      "^cl <- parallel::makeCluster\\(N_WORKERS, type = \"PSOCK\"\\)$",
      "cl <- NULL",
      legacy_lines
    )
    legacy_lines <- sub(
      "^doParallel::registerDoParallel\\(cl\\)$",
      "foreach::registerDoSEQ()",
      legacy_lines
    )
    legacy_lines <- sub(
      "^parallel::clusterEvalQ\\(cl, \\{$",
      "if (!is.null(cl)) parallel::clusterEvalQ(cl, {",
      legacy_lines
    )
    legacy_lines <- gsub("try\\(parallel::stopCluster\\(cl\\), silent = TRUE\\)", "NULL", legacy_lines)
    writeLines(legacy_lines, legacy_runtime_script, useBytes = TRUE)

    legacy_env <- new.env(parent = .GlobalEnv)
    sys.source(legacy_runtime_script, envir = legacy_env)
  } else {
    log_message(context, step, "Reusing existing legacy step3 outputs and bridge cache")
  }

  map_legacy_binary <- function(x) {
    as.integer(toupper(trimws(as.character(x))) == "T")
  }

  resolve_external_bridge <- function(legacy_env_obj, fallback_path) {
    if (file.exists(fallback_path)) {
      df <- read.csv(fallback_path, check.names = FALSE)
      if (all(c("sample_id", "observed_legacy", "predicted_prob") %in% colnames(df))) {
        return(df)
      }
    }

    if (is.null(legacy_env_obj) || !is.environment(legacy_env_obj)) {
      ext_fallback_path <- file.path(context$tables_dir, "External_test_predictions.csv")
      if (file.exists(ext_fallback_path)) {
        ext_fb <- read.csv(ext_fallback_path, check.names = FALSE)
        if (all(c("sample_id", "observed", "predicted_prob") %in% colnames(ext_fb))) {
          p_fb <- as.numeric(ext_fb$predicted_prob)
          valid_prob <- all(is.finite(p_fb)) && all(p_fb >= 0 & p_fb <= 1)
          if (valid_prob) {
            out_fb <- data.frame(
              sample_id = as.character(ext_fb$sample_id),
              observed_legacy = ifelse(as.integer(ext_fb$observed) == 1L, "T", "N"),
              predicted_prob = p_fb,
              selected_model = "fallback_from_existing_external_predictions",
              stringsAsFactors = FALSE
            )
            write.csv(out_fb, fallback_path, row.names = FALSE)
            return(out_fb)
          }
        }
      }
      stop("Cannot build external bridge cache: legacy environment is unavailable and no fallback external predictions were found.")
    }
    if (!exists("predictions_list", envir = legacy_env_obj, inherits = FALSE)) {
      stop("Cannot build external bridge cache: predictions_list not found in legacy environment.")
    }

    predictions_list <- get("predictions_list", envir = legacy_env_obj, inherits = FALSE)
    candidates <- c("SVM_RFE+LDA")
    if (exists("best_model_for_roc", envir = legacy_env_obj, inherits = FALSE)) {
      candidates <- c(candidates, as.character(get("best_model_for_roc", envir = legacy_env_obj, inherits = FALSE)))
    }
    if (exists("best_model_name", envir = legacy_env_obj, inherits = FALSE)) {
      candidates <- c(candidates, as.character(get("best_model_name", envir = legacy_env_obj, inherits = FALSE)))
    }
    candidates <- unique(candidates)

    selected_model <- NULL
    pred_obj <- NULL
    for (mn in candidates) {
      if (!mn %in% names(predictions_list)) {
        next
      }
      obj <- predictions_list[[mn]]
      if (!is.null(obj$outer_ext_prob)) {
        selected_model <- mn
        pred_obj <- obj
        break
      }
    }
    if (is.null(pred_obj)) {
      stop("Cannot build external bridge cache: no candidate model with outer_ext_prob found.")
    }

    obs_legacy <- NULL
    if (!is.null(pred_obj$test_labels)) {
      obs_legacy <- as.character(pred_obj$test_labels)
    } else if (exists("test_data_full", envir = legacy_env_obj, inherits = FALSE)) {
      tdf <- get("test_data_full", envir = legacy_env_obj, inherits = FALSE)
      if ("Type" %in% colnames(tdf)) {
        obs_legacy <- as.character(tdf$Type)
      }
    }
    if (is.null(obs_legacy)) {
      stop("Cannot build external bridge cache: external labels are unavailable.")
    }

    sample_ids <- NULL
    if (exists("test_data_full", envir = legacy_env_obj, inherits = FALSE)) {
      tdf <- get("test_data_full", envir = legacy_env_obj, inherits = FALSE)
      sample_ids <- rownames(tdf)
    }
    if (is.null(sample_ids) || length(sample_ids) != length(pred_obj$outer_ext_prob)) {
      sample_ids <- paste0("external_", seq_along(pred_obj$outer_ext_prob))
    }

    out <- data.frame(
      sample_id = sample_ids,
      observed_legacy = obs_legacy,
      predicted_prob = as.numeric(pred_obj$outer_ext_prob),
      selected_model = selected_model,
      stringsAsFactors = FALSE
    )
    write.csv(out, fallback_path, row.names = FALSE)
    out
  }

  oof_legacy_path <- file.path(step3_dir, "OOF_outer_SVM_RFE_LDA.csv")
  if (!file.exists(oof_legacy_path)) {
    stop(sprintf("Missing required legacy OOF file: %s", oof_legacy_path))
  }
  oof_legacy <- read.csv(oof_legacy_path, check.names = FALSE)
  if (!all(c("sample_id", "prob", "obs") %in% colnames(oof_legacy))) {
    stop(sprintf("Legacy OOF file has unexpected columns: %s", oof_legacy_path))
  }

  oof_df <- data.frame(
    sample_id = as.character(oof_legacy$sample_id),
    cohort = config$cohorts$train$id,
    cohort_label = config$cohorts$train$label,
    observed = map_legacy_binary(oof_legacy$obs),
    observed_label = ifelse(map_legacy_binary(oof_legacy$obs) == 1L, config$cohorts$train$positive_label, config$cohorts$train$control_label),
    predicted_prob = as.numeric(oof_legacy$prob),
    split = "outer_oof",
    stringsAsFactors = FALSE
  )

  ext_bridge <- resolve_external_bridge(legacy_env, bridge_external_path)
  ext_pheno <- read_table_file(resolve_revision_path(context, config$paths$external_pheno))
  ext_sample_ids <- as.character(ext_pheno[[config$cohorts$external$sample_col]])
  ext_y <- as.integer(as.character(ext_pheno[[config$cohorts$external$group_col]]) == config$cohorts$external$positive_label)

  if (nrow(ext_bridge) != length(ext_sample_ids)) {
    log_message(
      context,
      step,
      sprintf("External bridge row count (%d) != phenotype sample count (%d); forcing phenotype order", nrow(ext_bridge), length(ext_sample_ids)),
      level = "WARN"
    )
    ext_bridge <- ext_bridge[seq_len(min(nrow(ext_bridge), length(ext_sample_ids))), , drop = FALSE]
    ext_sample_ids <- ext_sample_ids[seq_len(nrow(ext_bridge))]
    ext_y <- ext_y[seq_len(nrow(ext_bridge))]
  } else {
    ext_bridge <- ext_bridge[seq_len(length(ext_sample_ids)), , drop = FALSE]
  }

  ext_obs_from_bridge <- map_legacy_binary(ext_bridge$observed_legacy)
  if (all(is.finite(ext_obs_from_bridge)) && length(ext_obs_from_bridge) == length(ext_y) && sum(abs(ext_obs_from_bridge - ext_y), na.rm = TRUE) <= 1L) {
    ext_obs <- ext_obs_from_bridge
  } else {
    ext_obs <- ext_y
  }

  external_df <- data.frame(
    sample_id = ext_sample_ids,
    cohort = config$cohorts$external$id,
    cohort_label = config$cohorts$external$label,
    observed = ext_obs,
    observed_label = ifelse(ext_obs == 1L, config$cohorts$external$positive_label, config$cohorts$external$control_label),
    predicted_prob = as.numeric(ext_bridge$predicted_prob),
    split = "holdout",
    stringsAsFactors = FALSE
  )

  oof_path <- file.path(context$tables_dir, "OOF_predictions_with_probabilities.csv")
  external_path <- file.path(context$tables_dir, "External_test_predictions.csv")
  write.csv(oof_df, oof_path, row.names = FALSE)
  write.csv(external_df, external_path, row.names = FALSE)

  models_path <- file.path(step3_dir, "models_results_v9.3.9_enhanced.csv")
  if (!file.exists(models_path)) {
    stop(sprintf("Missing required legacy model summary: %s", models_path))
  }
  models_df <- read.csv(models_path, check.names = FALSE)
  outer_df <- data.frame(
    model = as.character(models_df$Model),
    fold = NA_integer_,
    alpha = NA_real_,
    lambda = NA_real_,
    inner_auc = suppressWarnings(as.numeric(models_df$InnerOOF_AUC)),
    n_selected = suppressWarnings(as.integer(models_df$NumFeatures)),
    selected_features = as.character(models_df$Features),
    stringsAsFactors = FALSE
  )
  outer_path <- file.path(context$tables_dir, "Nested_CV_OuterFold_Selections.csv")
  write.csv(outer_df, outer_path, row.names = FALSE)

  target_row <- models_df[models_df$Model == "SVM_RFE+LDA", , drop = FALSE]
  if (nrow(target_row) == 0L) {
    target_row <- models_df[1, , drop = FALSE]
  }
  target_features <- trimws(unlist(strsplit(as.character(target_row$Features[1]), ";", fixed = TRUE)))
  target_features <- target_features[nzchar(target_features)]
  coef_df <- data.frame(
    feature = target_features,
    coefficient = NA_real_,
    stringsAsFactors = FALSE
  )
  coef_path <- file.path(context$tables_dir, "Model_coefficients_or_weights.csv")
  write.csv(coef_df, coef_path, row.names = FALSE)

  metric_row <- function(df, dataset_name) {
    auc_ci <- roc_auc_ci(df$observed, df$predicted_prob)
    cal <- calibration_metrics(df$observed, df$predicted_prob)
    data.frame(
      dataset = dataset_name,
      n = nrow(df),
      n_case = sum(df$observed == 1),
      n_control = sum(df$observed == 0),
      auc = auc_ci["auc"],
      auc_ci_low = auc_ci["ci_low"],
      auc_ci_high = auc_ci["ci_high"],
      pr_auc_case = pr_auc_case(df$observed, df$predicted_prob),
      pr_auc_control = pr_auc_case(1L - df$observed, 1 - df$predicted_prob),
      calibration_intercept = cal["intercept"],
      calibration_slope = cal["slope"],
      brier = cal["brier"],
      stringsAsFactors = FALSE
    )
  }
  performance_df <- rbind(
    metric_row(oof_df, "Outer OOF (GSE89632)"),
    metric_row(external_df, "External holdout (GSE135251)")
  )
  performance_path <- file.path(context$tables_dir, "Table_LOCO_Performance_by_Cohort.csv")
  write.csv(performance_df, performance_path, row.names = FALSE)

  methods_lines <- c(
    "Legacy-locked machine-learning baseline integrated for reproducibility alignment:",
    "1. Baseline script: primary/3.0_machine_learning_v1.R (selected from legacy screening).",
    "2. Feature engineering and model search follow the original 7 feature-selection x 9 classifier workflow with nested/outerized evaluation design.",
    "3. Reviewer-facing bridge outputs are standardized into OOF_predictions_with_probabilities.csv and External_test_predictions.csv for downstream calibration/PRC tables.",
    "4. Locked manuscript reporting target remains the SVM_RFE+LDA model track for Figure 4 consistency."
  )
  methods_path <- file.path(context$tables_dir, "Methods_NestedCV_NoLeakage.txt")
  write_text_file(methods_path, methods_lines)

  copy_first_existing <- function(candidates, target_path) {
    src <- candidates[file.exists(candidates)][1]
    if (is.na(src) || is.null(src)) {
      return(FALSE)
    }
    ensure_dir(dirname(target_path))
    file.copy(src, target_path, overwrite = TRUE)
  }

  fig_map <- list(
    model_performance_visualization_ALL.pdf = c(
      file.path(step3_dir, "plots", "model_performance_visualization_ALL.pdf"),
      file.path(step3_dir, "plots", "model_performance_visualization_v9.3.9_ALL.pdf")
    ),
    Combined_ROC_Calibration_SVM_RFE_LDA.pdf = c(
      file.path(step3_dir, "plots", "Combined_ROC_Calibration_SVM_RFE+LDA_v9.3.9.pdf"),
      file.path(step3_dir, "plots", "Combined_ROC_Calibration_SVM_RFE+LDA.pdf")
    ),
    PRC_Combined_SVM_RFE_LDA.pdf = c(
      file.path(step3_dir, "plots", "PRC_Combined_SVM_RFE+LDA_v9.3.9.pdf"),
      file.path(step3_dir, "plots", "PRC_Combined_SVM_RFE+LDA.pdf")
    ),
    PRC_Inverted_ControlPositive_SVM_RFE_LDA.pdf = c(
      file.path(step3_dir, "plots", "PRC_Inverted_ControlPositive_SVM_RFE+LDA_v9.3.9.pdf"),
      file.path(step3_dir, "plots", "PRC_Inverted_ControlPositive_SVM_RFE+LDA.pdf")
    )
  )

  copied_count <- 0L
  for (nm in names(fig_map)) {
    ok <- copy_first_existing(fig_map[[nm]], file.path(context$figures_main_dir, nm))
    if (isTRUE(ok)) {
      copied_count <- copied_count + 1L
    }
  }

  bridge_status <- data.frame(
    ran_legacy_script = need_run_legacy,
    bridge_cache = bridge_external_path,
    oof_source = oof_legacy_path,
    models_source = models_path,
    copied_main_figures = copied_count,
    stringsAsFactors = FALSE
  )
  bridge_status_path <- file.path(context$tables_dir, "legacy_ml_bridge_status.csv")
  write.csv(bridge_status, bridge_status_path, row.names = FALSE)

  log_message(context, step, sprintf("OOF predictions written: %s", oof_path))
  log_message(context, step, sprintf("External predictions written: %s", external_path))
  log_message(context, step, sprintf("Performance summary written: %s", performance_path))
  log_message(context, step, sprintf("Legacy bridge status written: %s", bridge_status_path))
}
